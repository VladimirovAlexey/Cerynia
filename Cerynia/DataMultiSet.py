"""
DataMultiSet class for the Cerynia library.

A sequence container of DataSets sharing the same process type.
Delegates all chi2 / matching / replica operations to constituent sets,
handling the slicing of a single flat theory vector back to per-set chunks.

Construction:
    multi = DataMultiSet([set1, set2, set3])
    multi = set1 + set2 + set3          # requires __add__ on DataSet

Access:
    multi[0]        # by index
    multi["CDF"]    # by name
    multi[0:2]      # slice → new DataMultiSet

Typical usage:
    multi = DataMultiSet([dy_set1, dy_set2])
    multi.prepare()
    theory = compute_theory(multi.df)   # flat vector, length = multi.numberOfPoints
    matched = multi.match(theory)
    chi2_total, chi2_list = multi.chi2(matched)
"""

import numpy as np
import pandas as pd

from .DataSet import DataSet


class DataMultiSet:

    def __init__(self, sets):
        if not sets:
            raise ValueError("DataMultiSet requires at least one DataSet")

        processTypes = {s.processType for s in sets}
        if len(processTypes) > 1:
            raise ValueError(
                f"All DataSets must share the same processType, got {processTypes}"
            )

        self._sets       = list(sets)
        self.processType = self._sets[0].processType

        # _i1[k], _i2[k]: start and end index of set k in a flat theory vector.
        # Set k occupies theory[_i1[k] : _i2[k]].
        self._i1 = []
        self._i2 = []
        i = 0
        for s in self._sets:
            self._i1.append(i)
            i += s.numberOfPoints
            self._i2.append(i)

        self.numberOfPoints = i

    # -- Representation --------------------------------------------------------

    def __repr__(self):
        state = "prepared" if all(s._prepared for s in self._sets) else "unprepared"
        return (f"<DataMultiSet: {len(self._sets)} sets, "
                f"{self.processType}, {self.numberOfPoints} points, {state}>")

    # -- Sequence interface ----------------------------------------------------

    def __len__(self):
        return len(self._sets)

    def __iter__(self):
        return iter(self._sets)

    def __getitem__(self, key):
        if isinstance(key, str):
            for s in self._sets:
                if s.name == key:
                    return s
            raise KeyError(f"No DataSet with name '{key}'")
        elif isinstance(key, slice):
            return DataMultiSet(self._sets[key])
        else:
            return self._sets[key]

    def __add__(self, other):
        if isinstance(other, DataSet):
            return DataMultiSet(self._sets + [other])
        elif isinstance(other, DataMultiSet):
            return DataMultiSet(self._sets + other._sets)
        return NotImplemented

    # -- Convenience -----------------------------------------------------------

    @property
    def df(self):
        """Concatenated DataFrame of all sets — useful for inspection and theory calls."""
        return pd.concat([s.df for s in self._sets], ignore_index=True)

    # -- Lifecycle -------------------------------------------------------------

    def prepare(self):
        """Call prepare() on all constituent DataSets. Returns self for chaining."""
        for s in self._sets:
            s.prepare()
        return self

    def cut(self, func):
        """
        Apply the same cut to all constituent DataSets.
        See DataSet.cut for the accepted forms of func.
        Sets left with zero points after the cut are dropped (with a printed
        notice); raises ValueError if every set ends up empty.
        Returns a new (unprepared) DataMultiSet.
        """
        kept = []
        for s in self._sets:
            cutSet = s.cut(func)
            if cutSet.numberOfPoints == 0:
                print(f"DataMultiSet.cut: all points removed from set '{s.name}' -- dropped.")
            else:
                kept.append(cutSet)
        if not kept:
            raise ValueError("DataMultiSet.cut: all points removed from every set")
        return DataMultiSet(kept)

    # -- Internal --------------------------------------------------------------

    def _slice(self, vector, k):
        """Extract the portion of a flat vector belonging to set k."""
        return vector[self._i1[k]:self._i2[k]]

    # -- Theory matching -------------------------------------------------------

    def match(self, theory):
        """
        Apply per-set matching (thFactor + normalization) to a flat theory vector.
        Returns a flat numpy array of the same length.
        """
        theory  = np.asarray(theory, dtype=float)
        matched = [self._sets[k].match(self._slice(theory, k))
                   for k in range(len(self._sets))]
        return np.concatenate(matched)

    # -- Chi2 and diagnostics --------------------------------------------------

    def chi2(self, theory):
        """
        Compute chi2 for all sets against a flat matched theory vector.
        Returns (chi2_total, [chi2_per_set]).
        """
        theory    = np.asarray(theory, dtype=float)
        chi2_list = [self._sets[k].chi2(self._slice(theory, k))
                     for k in range(len(self._sets))]
        return float(np.sum(chi2_list)), chi2_list

    def systematic_shift(self, theory):
        """
        Compute per-point systematic shifts for each set.
        Returns a list of arrays, one per set.
        """
        theory = np.asarray(theory, dtype=float)
        return [self._sets[k].systematic_shift(self._slice(theory, k))
                for k in range(len(self._sets))]

    def average_systematic_shift(self, theory):
        """
        Compute the mean systematic shift (as fraction of xSec) for each set.
        Returns a list of floats, one per set.
        """
        theory = np.asarray(theory, dtype=float)
        return [self._sets[k].average_systematic_shift(self._slice(theory, k))
                for k in range(len(self._sets))]

    def decompose_chi2(self, theory):
        """
        Decompose chi2 into correlated and uncorrelated parts for each set.
        Returns a list of (chi_D, chi_L, chi_total) tuples, one per set.
        """
        theory = np.asarray(theory, dtype=float)
        return [self._sets[k].decompose_chi2(self._slice(theory, k))
                for k in range(len(self._sets))]

    # -- Replica generation ----------------------------------------------------

    def generate_replica(self, include_norm_in_V=True, rng=None):
        """
        Generate a replica of all constituent DataSets.
        Returns a new prepared DataMultiSet.

        rng : numpy Generator (default: numpy.random.default_rng()).
              Pass numpy.random.default_rng(seed) for reproducibility, or share
              one generator across multiple calls to keep the random sequence coherent.
        """
        if rng is None:
            rng = np.random.default_rng()
        return DataMultiSet([s.generate_replica(include_norm_in_V, rng)
                             for s in self._sets])

    # -- Information -----------------------------------------------------------

    def info(self):
        """Print a summary of all constituent DataSets."""
        print(f"DataMultiSet : {self.processType}, "
              f"{len(self._sets)} sets, {self.numberOfPoints} points total")
        print(f"  {'Name':<30} {'N':>5}  {'Prepared':>8}")
        print("  " + "-" * 46)
        for s in self._sets:
            print(f"  {s.name:<30} {s.numberOfPoints:>5}  "
                  f"{'yes' if s._prepared else 'no':>8}")
