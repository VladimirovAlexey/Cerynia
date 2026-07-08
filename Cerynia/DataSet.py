"""
DataSet class for the Cerynia library.

Lifecycle:
    1. Loaded   — DataFrame + metadata, no matrices
    2. Cut      — filtered subset, still no matrices
    3. Prepared — covariance matrices computed, chi2 available

Typical usage:
    ds  = DataSet.from_csv("data/CDF.csv")
    sub = ds.cut(lambda df: df["qT_avg"] / df["Q_avg"] < 0.25)
    sub.prepare()
    chi2 = sub.chi2(matched_theory)
"""

import os
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from . import Point


class DataSet:

    def __init__(self, df, processType, name="", comment="", reference="",
                 normErr=None, isNormalized=False, normalizationMethod="integral"):
        if processType not in Point.PROCESS_TYPES:
            raise ValueError(f"processType must be one of {Point.PROCESS_TYPES}")

        self.df                  = df.reset_index(drop=True)
        self.processType         = processType
        self.name                = name
        self.comment             = comment
        self.reference           = reference
        self.normErr             = list(normErr) if normErr is not None else []
        self.isNormalized        = isNormalized
        self.normalizationMethod = normalizationMethod
        self._prepared           = False

        # Check for duplicate point ids
        if "id" in self.df.columns and self.df["id"].duplicated().any():
            dupes = self.df.loc[self.df["id"].duplicated(keep=False), "id"].unique().tolist()
            raise ValueError(f"DataSet '{name}' contains duplicate point ids: {dupes}")
            
    # -- Point-by-point construction -------------------------------------------

    @classmethod
    def empty(cls, processType, name="", comment="", reference="",
              normErr=None, isNormalized=False, normalizationMethod="integral"):
        """
        Create an empty DataSet with no points.
        Use add_point() to populate it, then prepare() when done.

        Example:
            ds = DataSet.empty("DY", name="CDF")
            for row in my_data:
                ds.add_point(row)
            ds.prepare()
        """
        return cls(
            df                  = pd.DataFrame(),
            processType         = processType,
            name                = name,
            comment             = comment,
            reference           = reference,
            normErr             = normErr,
            isNormalized        = isNormalized,
            normalizationMethod = normalizationMethod,
        )

    def add_point(self, point):
        """
        Validate and append a single point to the DataSet.
        point can be a dict or a pandas Series.
        Raises ValueError if a point with the same id already exists.
        Resets prepared state since the data has changed.
        """
        row_df = Point.validate(pd.DataFrame([point]), self.processType)

        point_id = row_df["id"].iloc[0]
        if len(self.df) > 0 and "id" in self.df.columns and point_id in self.df["id"].values:
            raise ValueError(
                f"A point with id '{point_id}' already exists in DataSet '{self.name}'"
            )

        self.df        = pd.concat([self.df, row_df], ignore_index=True)
        self._prepared = False


    # -- Representation --------------------------------------------------------

    def __repr__(self):
        state = "prepared" if self._prepared else "unprepared"
        return f"<DataSet: {self.name!r}, {self.processType}, {len(self.df)} points, {state}>"

    def __len__(self):
        return len(self.df)

    def __add__(self, other):
        # Lazy import to avoid circular dependency (DataMultiSet imports DataSet)
        from .DataMultiSet import DataMultiSet
        if isinstance(other, DataSet):
            return DataMultiSet([self, other])
        elif isinstance(other, DataMultiSet):
            return DataMultiSet([self] + list(other))
        return NotImplemented

    # -- Properties ------------------------------------------------------------

    @property
    def numberOfPoints(self):
        return len(self.df)

    @property
    def _uncorr_cols(self):
        # Sorted numerically (not alphabetically) so that uncorrErr_10 comes
        # after uncorrErr_9, preserving the intended order of error sources.
        cols = [c for c in self.df.columns if c.startswith("uncorrErr_")]
        return sorted(cols, key=lambda c: int(c.split("_")[1]))

    @property
    def _corr_cols(self):
        # Same numeric sort as _uncorr_cols — the order matters because C is
        # built column-by-column and must be consistent across prepare() calls.
        cols = [c for c in self.df.columns if c.startswith("corrErr_")]
        return sorted(cols, key=lambda c: int(c.split("_")[1]))

    def _require_prepared(self):
        if not self._prepared:
            raise RuntimeError(
                f"DataSet '{self.name}' is not prepared. Call prepare() first."
            )

    # -- Lifecycle -------------------------------------------------------------

    def prepare(self):
        """
        Compute and store various elements required to compute chi2,
        such as :   _variances = the list of uncorr^2 (N)
                    _C = the matrix of correlated uncertainties (N,nc)
                    _L = the Cholesky factor of the covariance matrix (N,N), 
                    _matrixA = matrix required to compute systematic shifts (nc,nc); see B4 in 1902.08474 
                    _matrixAinv = inverse of A.
        Returns self to allow chaining: sub.prepare().chi2(theory)
        """
        N = len(self.df)
        
        ## constructing the covariance matrix and related objects.
        ## the definition of V is taken from [1902.08474] sec.4, see also appendix A,B

        # Diagonal elements of the covariance matrix (sum of uncorrelated^2)        
        self._variances = self.df[self._uncorr_cols].pow(2).sum(axis=1).values  # (N,)
                
        # Build C: the (N x nc) matrix where row i is the full correlated error vector
        # for point i:  corr[i] = [corrErr_0[i], ..., xSec[i]*normErr_0, ...]
        # The covariance matrix is then V[i,j] = variances[i]*delta_ij + corr[i]·corr[j],
        # or in matrix form: V = diag(variances) + C @ C.T
        # Normalization uncertainties (stored as relative values in self.normErr) are
        # appended as xSec[i]*normErr_k — converting them to absolute per-point values
        # so they enter the formula on the same footing as corrErr columns.
        parts = []
        if self._corr_cols:
            parts.append(self.df[self._corr_cols].values)
        if self.normErr:
            parts.append(np.outer(self.df["xSec"].values, self.normErr))
        ### combines parts into the rows of C
        self._C = np.hstack(parts) if parts else np.zeros((N, 0))  # (N, nc)

        # Covariance matrix V[i,j] = variances[i]*delta_ij + corr[i]·corr[j]        
        # We store only its Cholesky factor L (V = L L.T), which is sufficient
        # for all chi2 computations and avoids keeping the full N×N matrix.
        V = np.diag(self._variances) + self._C @ self._C.T
        self._L = np.linalg.cholesky(V)

        # Matrix A = I + C_norm.T @ C_norm  (for systematic shift via nuisance parameters)
        # where C_norm[i]= corr[i]/sqrt(uncorr^2) (see B.4) in 1902.08474
        nc = self._C.shape[1]
        if nc > 0:
            C_norm = self._C / np.sqrt(self._variances)[:, None]
            self._matrixA    = np.eye(nc) + C_norm.T @ C_norm
            self._matrixAinv = np.linalg.inv(self._matrixA)
        else:
            self._matrixA    = np.eye(0)
            self._matrixAinv = np.eye(0)

        self._prepared = True
        return self

    def cut(self, func):
        """
        Return a new unprepared DataSet after filtering and optionally transforming points.

        func can be:
          - a boolean Series/array           : pure filter, selects Points where True
          - callable(Point) -> bool          : pure filter, applied Point by Point
          - callable(Point) -> (bool, Point) : filter + transform; the returned Point
                                               replaces the original (use this to modify
                                               point values as part of the selection)

        Examples:
            sub = ds.cut(ds.df["Q_avg"] > 10)
            sub = ds.cut(lambda p: p["xSec"] > 0)
            sub = ds.cut(lambda p: (p["Q_avg"] > 10, p))
            def my_cut(p):
                keep = p["xSec"] / p["uncorrErr_0"] > 2   # signal-to-noise cut
                pN = p.copy()
                pN["thFactor"] *= 0.5                         # rescale on the fly
                return keep, pN
            sub = ds.cut(my_cut)
        """
        if not callable(func):
            # Plain boolean mask — simple filter, no transformation
            new_df = self.df[func]
        else:
            # Apply func to each row; detect whether it returns (bool, row) or just bool
            results = [func(row) for _, row in self.df.iterrows()]
            if isinstance(results[0], tuple):
                # (bool, modified_row) interface
                rows = [row for include, row in results if include]
            else:
                # bool-only interface — keep original row unchanged
                rows = [row for include, (_, row) in zip(results, self.df.iterrows()) if include]
            new_df = pd.DataFrame(rows)

        return DataSet(
            df                  = new_df,
            processType         = self.processType,
            name                = self.name,
            comment             = self.comment,
            reference           = self.reference,
            normErr             = self.normErr.copy(),
            isNormalized        = self.isNormalized,
            normalizationMethod = self.normalizationMethod,
        )

    # -- Theory matching -------------------------------------------------------

    def match(self, theory):
        """
        Apply thFactor and optional normalization to a raw theory vector.
        Returns a numpy array aligned with self.df.

        For bestChi2 normalization, prepare() must have been called.
        """
        matched = np.asarray(theory, dtype=float) * self.df["thFactor"].values

        if not self.isNormalized:
            return matched

        if self.normalizationMethod == "integral":
            if self.processType == "DY":
                width = self.df["qT_max"].values - self.df["qT_min"].values
            elif self.processType == "SIDIS":
                width = self.df["pT_max"].values - self.df["pT_min"].values
            else:
                width = np.ones(len(self.df))
            norm_exp = np.dot(width, self.df["xSec"].values)
            norm_th  = np.dot(width, matched)
            return matched * (norm_exp / norm_th)

        if self.normalizationMethod == "bestChi2":
            self._require_prepared()
            return matched * self._best_norm(matched)

        raise ValueError(f"Unknown normalizationMethod: {self.normalizationMethod!r}")

    # -- Chi2 and diagnostics --------------------------------------------------

    def chi2(self, theory):
        """
        Compute chi2 = (theory - xSec)^T V^{-1} (theory - xSec).
        theory should already be matched via match().
        """
        self._require_prepared()
        x  = np.asarray(theory, dtype=float) - self.df["xSec"].values
        Lx = solve_triangular(self._L, x, lower=True)
        return float(np.dot(Lx, Lx))

    def systematic_shift(self, theory):
        """
        Compute per-point systematic shifts via nuisance parameters.
        Returns a numpy array of shifts aligned with self.df.
        See eqn, B3-B5 in 1902.08474
        """
        self._require_prepared()
        if self._C.shape[1] == 0:
            return np.zeros(len(self.df))
        diff  = self.df["xSec"].values - np.asarray(theory, dtype=float)
        rho   = (diff / self._variances) @ self._C
        lambd = self._matrixAinv @ rho
        return self._C @ lambd

    def average_systematic_shift(self, theory):
        """Return the mean systematic shift as a fraction of xSec (e.g. 0.03 = 3%)."""
        shift = self.systematic_shift(theory)
        xSec  = self.df["xSec"].values
        mask  = xSec != 0
        return float(np.mean(shift[mask] / xSec[mask]))

    def decompose_chi2(self, theory):
        """
        Decompose chi2 into correlated and uncorrelated parts.
        Returns (chi_D, chi_L, chi_total).
        """
        self._require_prepared()
        theory = np.asarray(theory, dtype=float)
        shift  = self.systematic_shift(theory)

        chi_D = float(np.sum(
            (self.df["xSec"].values - theory - shift) ** 2 / self._variances
        ))

        if self._C.shape[1] > 0:
            diff  = self.df["xSec"].values - theory
            rho   = (diff / self._variances) @ self._C
            lambd = self._matrixAinv @ rho
            chi_L = float(np.dot(lambd, lambd))
        else:
            chi_L = 0.0

        return chi_D, chi_L, chi_D + chi_L

    def _best_norm(self, theory):
        """Normalization n minimizing chi2: n = (xSec · V⁻¹ · t) / (t · V⁻¹ · t)."""
        Ls = solve_triangular(self._L, self.df["xSec"].values, lower=True)
        Lt = solve_triangular(self._L, theory,                 lower=True)
        return float(np.dot(Ls, Lt) / np.dot(Lt, Lt))

    # -- Replica generation ----------------------------------------------------

    def generate_replica(self, include_norm_in_V=True, rng=None):
        """
        Generate a pseudo-data replica by Gaussian fluctuation of the data.
        Recipe from arXiv:0808.1231 sec. 2.4.
        Note, that all normalization uncertanties are considered as "absolute" according to sec. 2.2

        include_norm_in_V : include normErr in the replica's covariance matrix
        rng               : numpy Generator (default: numpy.random.default_rng()).
                            Pass numpy.random.default_rng(seed) for reproducibility,
                            or share one generator across multiple calls to correlate
                            the random sequences.
        """
        if rng is None:
            rng = np.random.default_rng()

        df = self.df.copy()

        # Common normalization fluctuation: prod(1 + N(0,1) * normErr_i)
        norm_factor = (
            float(np.prod(1 + rng.standard_normal(len(self.normErr)) * np.array(self.normErr)))
            if self.normErr else 1.0
        )

        # Fluctuate xSec
        xSec = df["xSec"].values.copy().astype(float)

        for col in self._uncorr_cols:                              # uncorrelated: one RND per point
            xSec += rng.standard_normal(len(df)) * df[col].values

        if self._corr_cols:                                        # correlated: one RND per source
            corr_rnd = rng.standard_normal(len(self._corr_cols))
            xSec += df[self._corr_cols].values @ corr_rnd

        xSec *= norm_factor
        df["xSec"] = xSec

        # Rescale all errors by the normalization factor
        for col in self._uncorr_cols + self._corr_cols:
            df[col] = df[col].values * norm_factor

        return DataSet(
            df                  = df,
            processType         = self.processType,
            name                = self.name + "(rep)",
            comment             = self.comment,
            reference           = self.reference,
            normErr             = self.normErr if include_norm_in_V else [],
            isNormalized        = self.isNormalized,
            normalizationMethod = self.normalizationMethod,
        ).prepare()

    # -- Information -----------------------------------------------------------

    def info(self):
        """Print a human-readable summary of the dataset."""
        print(f"DataSet      : {self.name}")
        if self.reference:
            print(f"Reference    : {self.reference}")
        if self.comment:
            print(f"Comment      : {self.comment}")
        print(f"Process type : {self.processType}")
        print(f"Points       : {len(self.df)}")
        print(f"Uncorr. err  : {len(self._uncorr_cols)}")
        print(f"Corr. err    : {len(self._corr_cols)}")
        print(f"Norm. err    : {len(self.normErr)}")
        norm_str = f"yes ({self.normalizationMethod})" if self.isNormalized else "no"
        print(f"Normalized   : {norm_str}")
        print(f"Prepared     : {self._prepared}")

    # -- I/O -------------------------------------------------------------------

    @classmethod
    def from_csv(cls, path):
        """
        Load a DataSet from a CSV file.
        Metadata is read from leading comment lines of the form:
            # key: value
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"Data file not found: {path}")

        meta = {}
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                key, _, val = line[2:].partition(":")
                meta[key.strip()] = val.strip()

        normErr     = [float(x) for x in meta.get("normErr", "").split(",") if x.strip()]
        processType = meta.get("processType", "")

        df = pd.read_csv(path, comment="#")
        df = Point.validate(df, processType)

        return cls(
            df                  = df,
            processType         = processType,
            name                = meta.get("name", ""),
            comment             = meta.get("comment", ""),
            reference           = meta.get("reference", ""),
            normErr             = normErr,
            isNormalized        = meta.get("isNormalized", "False") == "True",
            normalizationMethod = meta.get("normalizationMethod", "integral"),
        )

    def save_csv(self, path):
        """
        Save the DataSet to a CSV file.
        Metadata is written as leading comment lines; the DataFrame follows.
        """
        with open(path, "w") as f:
            f.write(f"# name: {self.name}\n")
            f.write(f"# comment: {self.comment}\n")
            f.write(f"# reference: {self.reference}\n")
            f.write(f"# processType: {self.processType}\n")
            f.write(f"# isNormalized: {self.isNormalized}\n")
            f.write(f"# normalizationMethod: {self.normalizationMethod}\n")
            f.write(f"# normErr: {','.join(str(e) for e in self.normErr)}\n")
        self.df.to_csv(path, mode="a", index=False)
