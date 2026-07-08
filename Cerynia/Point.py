"""
Point schema definitions and validation for the Cerynia library.

A "point" is a row in a DataSet's DataFrame. This module defines:
  - PROCESS_TYPES : valid process type identifiers
  - REQUIRED      : required columns per process type
  - validate()    : fill defaults and check required columns
  - create()      : build a validated single-row DataFrame
  - schema()      : print a human-readable schema summary
"""

import pandas as pd

PROCESS_TYPES = ("DY", "SIDIS", "G2", "D2")

# -- Required columns -----------------------------------------------------------

_COMMON_REQUIRED = [
    "id", "xSec", "s", "Q_min", "Q_max", "thFactor",
    "ps_def", "h_1", "proc_id",
]

_EXTRA_REQUIRED = {
    "DY":    ["h_2", "y_min", "y_max", "qT_min", "qT_max", "includeCuts"],
    "SIDIS": ["h_2", "x_min", "x_max", "z_min", "z_max", "pT_min", "pT_max", "includeCuts"],
    "G2":    ["x_min", "x_max"],
    "D2":    [],
}

REQUIRED = {pt: _COMMON_REQUIRED + _EXTRA_REQUIRED[pt] for pt in PROCESS_TYPES}

# -- Optional columns with defaults --------------------------------------------
# Ordered list of (column, value_or_callable).
# Callables receive the current DataFrame and return a Series.
# Order matters: later entries may depend on earlier computed columns.

_DEFAULTS = {
    "DY": [
        ("Q_avg",       lambda df: (df["Q_min"] + df["Q_max"]) / 2),
        ("y_avg",       lambda df: (df["y_min"] + df["y_max"]) / 2),
        ("qT_avg",      lambda df: (df["qT_min"] + df["qT_max"]) / 2),
        ("cutParams_0", 0.0),
        ("cutParams_1", 0.0),
        ("cutParams_2", -100.0),
        ("cutParams_3", 100.0),
    ],
    "SIDIS": [
        ("Q_avg",       lambda df: (df["Q_min"] + df["Q_max"]) / 2),
        ("x_avg",       lambda df: (df["x_min"] + df["x_max"]) / 2),
        ("z_avg",       lambda df: (df["z_min"] + df["z_max"]) / 2),
        ("pT_avg",      lambda df: (df["pT_min"] + df["pT_max"]) / 2),        
        ("cutParams_0", 0.0),
        ("cutParams_1", 1.0),
        ("cutParams_2", 0.0),
        ("cutParams_3", 100000.0),
        ("M_target",    0.938),
        ("M_product",   0.139),
    ],
    "G2": [
        ("Q_avg",  lambda df: (df["Q_min"] + df["Q_max"]) / 2),
        ("x_avg",  lambda df: (df["x_min"] + df["x_max"]) / 2),
    ],
    "D2": [
        ("Q_avg",  lambda df: (df["Q_min"] + df["Q_max"]) / 2),
    ],
}

# Optional weight-process columns (DY/SIDIS only, present only when needed)
_WEIGHT_COLUMNS = {
    "DY":    ["ps_def_weight", "h_1_weight", "h_2_weight", "proc_id_weight"],
    "SIDIS": ["ps_def_weight", "h_1_weight", "h_2_weight", "proc_id_weight"],
    "G2":    [],
    "D2":    [],
}

# -- Internal helpers -----------------------------------------------------------

def _fill_defaults(df, processType):
    df = df.copy()
    for col, default in _DEFAULTS[processType]:
        if col not in df.columns:
            df[col] = default(df) if callable(default) else default
    return df


def _check_bin_order(df, col_min, col_max):
    bad = df[col_min] > df[col_max]
    if bad.any():
        bad_ids = df.loc[bad, "id"].tolist()
        raise ValueError(
            f"Bin ordering violated: '{col_min}' > '{col_max}' for points: {bad_ids}"
        )


# -- Public API -----------------------------------------------------------------

def validate(df, processType):
    """
    Validate a DataFrame of points against the schema for processType.

    Fills missing optional columns with defaults, then checks required columns
    and bin orderings. Returns the validated (possibly augmented) DataFrame.
    Raises ValueError on any violation.
    """
    if processType not in PROCESS_TYPES:
        raise ValueError(f"processType must be one of {PROCESS_TYPES}, got '{processType}'")

    missing = [c for c in REQUIRED[processType] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for {processType}: {missing}")

    df = _fill_defaults(df, processType)

    _check_bin_order(df, "Q_min", "Q_max")
    if processType == "DY":
        _check_bin_order(df, "y_min",  "y_max")
        _check_bin_order(df, "qT_min", "qT_max")
    elif processType == "SIDIS":
        _check_bin_order(df, "x_min",  "x_max")
        _check_bin_order(df, "z_min",  "z_max")
        _check_bin_order(df, "pT_min", "pT_max")
    elif processType == "G2":
        _check_bin_order(df, "x_min",  "x_max")

    return df


def create(processType, **kwargs):
    """
    Build a validated single-row DataFrame for one data point.

    Example
    -------
    p = Point.create("DY",
            id="CDF-1", ps_def=1, h_1=2212, h_2=-2212, proc_id=1,
            xSec=1.23, s=3600., Q_min=66., Q_max=116., thFactor=1.,
            y_min=-1., y_max=1., qT_min=0., qT_max=5.,
            includeCuts=False,
            uncorrErr_0=0.05)
    """
    return validate(pd.DataFrame([kwargs]), processType)


def schema(processType):
    """Print a human-readable summary of columns for processType."""
    if processType not in PROCESS_TYPES:
        raise ValueError(f"processType must be one of {PROCESS_TYPES}")

    print(f"Process type : {processType}")
    print("  Required columns:")
    for c in REQUIRED[processType]:
        print(f"    {c}")
    print("  Optional columns (with defaults):")
    for col, default in _DEFAULTS[processType]:
        label = "computed" if callable(default) else str(default)
        print(f"    {col:<20} [{label}]")
    if _WEIGHT_COLUMNS[processType]:
        print("  Optional weight-process columns:")
        for c in _WEIGHT_COLUMNS[processType]:
            print(f"    {c}")
