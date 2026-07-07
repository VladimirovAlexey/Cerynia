"""
harpy/artemide interface for the Cerynia library.

Handles all four process types by dispatching on data.processType:
    DY    — harpy.DY.xSecList   (bin-integrated, optional partitioning)
    SIDIS — harpy.SIDIS.xSecList (bin-integrated, optional partitioning)
    G2    — harpy.G2List         (evaluated at bin-center averages)
    D2    — harpy.D2List         (evaluated at bin-center averages)

Public API:
    xsec(data, method)              — raw theory vector from harpy
    chi2(data, method)              — xsec + match + chi2
    print_chi2_table(data, ...)     — formatted chi2/N table
    print_per_point_chi2(data, ...) — per-point diagonal chi2
"""

import time
import numpy as np
import harpy

from .DataSet      import DataSet
from .DataMultiSet import DataMultiSet


# ── Column extraction helpers ─────────────────────────────────────────────────

def _process_list(df, processType):
    """Process code list per point.
    DY/SIDIS: [ps_def, h_1, h_2, proc_id]  (4 integers)
    G2/D2:    [ps_def, h_1, proc_id]        (3 integers)
    """
    if processType in ("DY", "SIDIS"):
        return df[["ps_def", "h_1", "h_2", "proc_id"]].values.tolist()
    else:
        return df[["ps_def", "h_1", "proc_id"]].values.tolist()

def _cutparams_list(df):
    return df[["cutParams_0", "cutParams_1", "cutParams_2", "cutParams_3"]].values.tolist()

def _masses_list(df):
    return df[["M_target", "M_product"]].values.tolist()

def _zero_size_bins(df, col):
    """Zero-size bin [[v, v], ...] for 'central' method evaluation at bin centers."""
    return [[v, v] for v in df[col]]


# ── Core cross-section computation ────────────────────────────────────────────

def _xsec_df(df, processType, method):
    """
    Compute raw cross-sections for a flat DataFrame of points.
    Returns a numpy array of length len(df).
    """

    # G2 and D2: always evaluated at bin-center averages, no bin-integration variants.
    if processType == "G2":
        if method != "default":
            raise ValueError(f"G2 only supports method='default', got '{method}'")
        return np.array(harpy.G2List(
            df["x_avg"].tolist(),
            df["Q_avg"].tolist(),
            _process_list(df, processType),
        ))

    if processType == "D2":
        if method != "default":
            raise ValueError(f"D2 only supports method='default', got '{method}'")
        return np.array(harpy.D2List(
            df["Q_avg"].tolist(),
            _process_list(df, processType),
        ))

    # DY and SIDIS: full bin-integrated methods.
    if method in ("default", "noPartitioning"):
        partitioning = (method == "default")

        if processType == "DY":
            return np.array(harpy.DY.xSecList(
                _process_list(df, processType),
                df["s"].tolist(),
                df[["qT_min", "qT_max"]].values.tolist(),
                df[["Q_min",  "Q_max" ]].values.tolist(),
                df[["y_min",  "y_max" ]].values.tolist(),
                df["includeCuts"].tolist(),
                _cutparams_list(df),
                partitioning,
            ))

        if processType == "SIDIS":
            return np.array(harpy.SIDIS.xSecList(
                _process_list(df, processType),
                df["s"].tolist(),
                df[["pT_min", "pT_max"]].values.tolist(),
                df[["z_min",  "z_max" ]].values.tolist(),
                df[["x_min",  "x_max" ]].values.tolist(),
                df[["Q_min",  "Q_max" ]].values.tolist(),
                df["includeCuts"].tolist(),
                _cutparams_list(df),
                _masses_list(df),
                partitioning,
            ))

    if method == "central":
        # Evaluate at bin centers using zero-size bins [avg, avg]
        # (current artemide convention for point-like integration)
        if processType == "DY":
            return np.array(harpy.DY.xSecList(
                _process_list(df, processType),
                df["s"].tolist(),
                _zero_size_bins(df, "qT_avg"),
                _zero_size_bins(df, "Q_avg"),
                _zero_size_bins(df, "y_avg"),
                df["includeCuts"].tolist(),
                _cutparams_list(df),
            ))

        if processType == "SIDIS":
            return np.array(harpy.SIDIS.xSecList(
                _process_list(df, processType),
                df["s"].tolist(),
                _zero_size_bins(df, "pT_avg"),
                _zero_size_bins(df, "z_avg"),
                _zero_size_bins(df, "x_avg"),
                _zero_size_bins(df, "Q_avg"),
                [False] * len(df),
                _cutparams_list(df),
                _masses_list(df),
                False,
            ))

    raise ValueError(
        f"processType='{processType}', method='{method}': "
        "valid methods are 'default', 'noPartitioning', 'central' (DY/SIDIS) "
        "or 'default' (G2/D2)."
    )


# ── Public API ────────────────────────────────────────────────────────────────

def xsec(data, method="default"):
    """
    Compute raw cross-section / observable values from harpy.
    Returns a flat numpy array of length data.numberOfPoints.

    This is the RAW theory output — call data.match(result) afterwards to apply
    thFactor and normalization before passing to chi2.

    method (DY / SIDIS):
        'default'        — full bin integration with qT partitioning
        'noPartitioning' — full bin integration without partitioning
        'central'        — evaluate at bin centers only (zero-size bins)
    method (G2 / D2):
        'default'        — evaluate at bin-center averages (x_avg, Q_avg)
    """
    if not isinstance(data, (DataSet, DataMultiSet)):
        raise TypeError("data must be a DataSet or DataMultiSet")
    return _xsec_df(data.df, data.processType, method)


def chi2(data, method="default"):
    """
    Compute raw cross-sections, match to data, and evaluate chi2.
    Returns (chi2_total, [chi2_per_set]).
    data must be prepared.
    """
    XX = xsec(data, method)
    YY = data.match(XX)

    if isinstance(data, DataSet):
        result = data.chi2(YY)
        return result, [result]
    return data.chi2(YY)


# ── Diagnostic printing ───────────────────────────────────────────────────────

def print_chi2_table(data, method="default", decompose=False, sys_shift=True):
    """
    Compute and print a chi2/N summary table.

    method    : xsec computation method (see xsec())
    decompose : show chi_D^2/N (uncorrelated) and chi_L^2/N (correlated) separately
    sys_shift : show average systematic shift per set (in %)
    """
    t0 = time.time()

    XX = xsec(data, method)
    YY = data.match(XX)

    if isinstance(data, DataSet):
        sets       = [data]
        chi2_list  = [data.chi2(YY)]
        dec_list   = [data.decompose_chi2(YY)]            if decompose  else None
        shift_list = [data.average_systematic_shift(YY)]  if sys_shift  else None
    else:
        sets       = list(data)
        _, chi2_list = data.chi2(YY)
        dec_list   = data.decompose_chi2(YY)          if decompose  else None
        shift_list = data.average_systematic_shift(YY) if sys_shift  else None

    dt = time.time() - t0
    w  = max(len(s.name) for s in sets)

    header = f"{'name':{w}} | {'N':>5} |"
    sep    = f"{'':-<{w}}-+-{'':->5}-+-"
    if decompose:
        header += f" {'chiD^2/N':>10} | {'chiL^2/N':>10} | {'chi^2/N':>10} |"
        sep    += f"-{'':->10}-+-{'':->10}-+-{'':->10}-+-"
    else:
        header += f" {'chi^2/N':>10} |"
        sep    += f"-{'':->10}-+-"
    if sys_shift:
        header += f" {'sys.shift%':>10} |"
        sep    += f"-{'':->10}-+-"

    print(header)
    print(sep)

    for i, s in enumerate(sets):
        N    = max(s.numberOfPoints, 1)
        line = f"{s.name:{w}} | {s.numberOfPoints:>5} |"
        if decompose:
            d = dec_list[i]
            line += f" {d[0]/N:>10.3f} | {d[1]/N:>10.3f} | {d[2]/N:>10.3f} |"
        else:
            line += f" {chi2_list[i]/N:>10.3f} |"
        if sys_shift:
            line += f" {shift_list[i]*100:>10.3f} |"
        print(line)

    if len(sets) > 1:
        N_total    = max(data.numberOfPoints, 1)
        chi2_total = sum(chi2_list)
        print(sep)
        line = f"{'Total':{w}} | {data.numberOfPoints:>5} |"
        if decompose:
            totals = np.sum(dec_list, axis=0)
            line += (f" {totals[0]/N_total:>10.3f} |"
                     f" {totals[1]/N_total:>10.3f} |"
                     f" {totals[2]/N_total:>10.3f} |")
        else:
            line += f" {chi2_total/N_total:>10.3f} |"
        print(line)

    print(f"Computation time: {dt:.2f} s")


def print_per_point_chi2(data, method="default", min_chi2=0.):
    """
    Print the diagonal chi2 contribution per point (uncorrelated errors only).
    Points with chi2 < min_chi2 are suppressed.
    """
    XX = xsec(data, method)
    YY = data.match(XX)

    if isinstance(data, DataSet):
        sets_and_theory = [(data, YY)]
    else:
        sets_and_theory = [(data[k], data._slice(YY, k)) for k in range(len(data))]

    for ds, th in sets_and_theory:
        if isinstance(data, DataMultiSet):
            print(f"\n── {ds.name} ──")

        variances = ds.df[ds._uncorr_cols].pow(2).sum(axis=1).values
        dchi2     = (ds.df["xSec"].values - th) ** 2 / variances

        w = max(ds.df["id"].str.len().max(), 10)
        print(f"  {'id':{w}} | {'dChi^2':>10}")
        print(f"  {'':-<{w}}-+-{'':-<10}")

        shown = 0
        for i in range(len(ds.df)):
            if dchi2[i] >= min_chi2:
                print(f"  {ds.df.iloc[i]['id']:{w}} | {dchi2[i]:>10.3f}")
                shown += 1
        if shown == 0:
            print(f"  All points have dChi^2 < {min_chi2}")
