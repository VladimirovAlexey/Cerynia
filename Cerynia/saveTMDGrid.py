"""
saveTMDGrid.py -- write TMD/TMDFF grids (as returned by harpy/artemide) to the
text grid format used by TMDlib and related codes. The file format itself is a
field standard and must not change:

    xg: [...]          (or Qg:/bg:/qToQg: as applicable, one line per axis)
    TMDs: {-5: [...], -4: [...], ..., 5: [...]}

where each leaf is 'x * TMD(x, ...)' formatted with '{:g}'.format, and the
dict keys are PDG-style flavor codes (b, c, s, u, d, [g], d, u, s, c, b).

Public API:
    save_grid_optimal(path, ...)    -- b-space, optimal (or single-Q) grid
    save_grid_optimal_kT(path, ...) -- kT-space, optimal (or single-Q) grid
    save_grid_Q(path, ...)          -- b-space grid over a Q range
    save_grid_kT(path, ...)         -- kT-space grid over a Q range
"""

import time
import numpy as np
import harpy


# -- Default grids ---------------------------------------------------------

### Q[GeV]
Q_RANGE_DEFAULT = [1., 1.11803, 1.22474, 1.4, 1.58114, 1.78885, 2., 2.23607, 2.52982, 2.82843,
    3.16228, 3.4641, 4.75, 5.09902, 6.32456, 7.1, 8., 10., 11.1803, 12.2475,
    14., 15.8114, 17.8885, 20., 22.3607, 25.2982, 28.2843, 31.6228, 34.641, 47.5,
    50.9902, 63.2456, 71, 80, 100, 111.803, 122.475, 140, 158.114, 178.885,
    200.]
### x for TMDPDFs
X_RANGE_PDF_DEFAULT = [0.00001, 0.00002, 0.00004, 0.00006, 0.00008, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008,
    0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055,
    0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.00925, 0.0095, 0.00975,
    0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055,
    0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.0925, 0.095, 0.0975,
    0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
    0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975,
    1]
### x for TMDFFs
X_RANGE_FF_DEFAULT = [0.05, 0.055, 0.06, 0.065, 0.07, 0.08, 0.09, 0.1,
    0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.]
### R = qT/Q for grids in the momentum space
R_RANGE_DEFAULT = [0.0001, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.04, 0.05,
    0.06, 0.07, 0.08, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225,
    0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
    0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
    1.7, 1.8, 1.9, 2.001]
B_RANGE_DEFAULT = [0., 0.01, 0.025, 0.05, 0.1,
    0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
    0.6, 0.65, 0.7, 0.75, 0.8, 0.86, 0.93, 1.,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.65, 1.8, 1.95, 2.1,
    2.3, 2.5, 2.75, 3., 3.25, 3.5, 4., 4.5, 5., 6, 8, 10.]
KT_RANGE_DEFAULT = [0.00001, 0.01, 0.025, 0.05, 0.1,
    0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
    0.6, 0.65, 0.7, 0.75, 0.8, 0.86, 0.93, 1.,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
    2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9,
    4.1, 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5,
    10., 11., 12., 13., 14., 15., 20., 25.]


# -- harpy dispatch ----------------------------------------------------------

_GETTERS = {
    "uTMDPDF":      "get_uTMDPDF",
    "uTMDFF":       "get_uTMDFF",
    "SiversTMDPDF": "get_SiversTMDPDF",
    "wgtTMDPDF":    "get_wgtTMDPDF",
}

_FLAVORS      = tuple(range(-5, 6))          # -5..5, includes gluon (0)
_FLAVORS_NOG  = tuple(f for f in _FLAVORS if f != 0)


def _getter(PDF, kT):
    try:
        name = _GETTERS[PDF]
    except KeyError:
        raise ValueError(f"wrong type for PDF: '{PDF}'; must be one of {tuple(_GETTERS)}")
    return getattr(harpy, name + "_kT" if kT else name)


def _default_Xrange(Xrange, PDF):
    """Xrange=None picks X_RANGE_FF_DEFAULT for uTMDFF, X_RANGE_PDF_DEFAULT otherwise."""
    if Xrange is not None:
        return Xrange
    return X_RANGE_FF_DEFAULT if PDF == "uTMDFF" else X_RANGE_PDF_DEFAULT


# -- grid computation ---------------------------------------------------------

def _build_grid(ranges, TMDval_of, flavors, fill_flavors):
    """
    ranges       : list of grid arrays, one per axis, outermost first
    TMDval_of    : callable(*axis_values) -> 11-element list of x*TMD (flavor -5..5)
    flavors      : flavor keys present in the output dict
    fill_flavors : subset of `flavors` actually computed (others stay 0.)
    """
    shape  = tuple(len(r) for r in ranges)
    values = {f: np.full(shape, 0., dtype=object) for f in flavors}

    for idx in np.ndindex(shape):
        axis_vals = tuple(float(ranges[a][idx[a]]) for a in range(len(ranges)))
        TMDval    = TMDval_of(*axis_vals)
        for f in fill_flavors:
            values[f][idx] = '{:g}'.format(TMDval[f + 5])

    return {f: v.tolist() for f, v in values.items()}


def _write(path, header, valuesList, t0):
    with open(path, 'w') as outfile:
        for key, grid in header:
            outfile.write(f"{key}: {grid}\n")
        outfile.write("TMDs: " + str(valuesList).replace("'", "") + "\n")
    print('Grid written at  : ', path)
    print('Computation time  : ', time.time() - t0, ' sec.')


# -- Public API ----------------------------------------------------------------

def save_grid_optimal(path, Xrange=None, Brange=B_RANGE_DEFAULT, PDF="uTMDPDF", h=1, Q=-1.):
    """
    Save the optimal (or single-Q, if Q>0) TMDPDF/TMDFF grid in b-space, as
    returned by harpy (artemide must already be set up).

    PDF    : which TMD to save ('uTMDPDF', 'uTMDFF', 'SiversTMDPDF', 'wgtTMDPDF')
    h      : hadron number
    Q      : scale to evaluate at; if negative, harpy's optimal grid is used
    Xrange : if omitted, defaults to X_RANGE_FF_DEFAULT for PDF='uTMDFF',
             X_RANGE_PDF_DEFAULT otherwise
    """
    t0     = time.time()
    Xrange = _default_Xrange(Xrange, PDF)
    getter = _getter(PDF, kT=False)

    def TMDval_of(xval, bval):
        if xval > 0.999:
            return [0.] * 11
        return [xval * v for v in getter(xval, bval, h, includeGluon=False, mu=Q)]

    valuesList = _build_grid([Xrange, Brange], TMDval_of, _FLAVORS_NOG, _FLAVORS_NOG)
    _write(path, [("xg", Xrange), ("bg", Brange)], valuesList, t0)


def save_grid_optimal_kT(path, Xrange=None, KTrange=KT_RANGE_DEFAULT, PDF="uTMDPDF", h=1, Q=-1.):
    """
    Save the optimal (or single-Q, if Q>0) TMDPDF/TMDFF grid in kT-space, as
    returned by harpy (artemide must already be set up).

    PDF    : which TMD to save ('uTMDPDF', 'uTMDFF', 'SiversTMDPDF', 'wgtTMDPDF')
    h      : hadron number
    Q      : scale to evaluate at; if negative, harpy's optimal grid is used
    Xrange : if omitted, defaults to X_RANGE_FF_DEFAULT for PDF='uTMDFF',
             X_RANGE_PDF_DEFAULT otherwise
    """
    t0     = time.time()
    Xrange = _default_Xrange(Xrange, PDF)
    getter = _getter(PDF, kT=True)

    def TMDval_of(xval, ktval):
        if xval > 0.999:
            return [0.] * 11
        return [xval * v for v in getter(xval, ktval, h, includeGluon=False, mu=Q)]

    valuesList = _build_grid([Xrange, KTrange], TMDval_of, _FLAVORS_NOG, _FLAVORS_NOG)
    _write(path, [("xg", Xrange), ("bg", KTrange)], valuesList, t0)


def save_grid_Q(path, Qrange=Q_RANGE_DEFAULT, Xrange=None, Brange=B_RANGE_DEFAULT,
                PDF="uTMDPDF", h=1, includeGluon=False):
    """
    Save the TMDPDF/TMDFF grid in b-space over a range of Q, Q^2, as returned
    by harpy (artemide must already be set up).

    PDF          : which TMD to save ('uTMDPDF', 'uTMDFF', 'SiversTMDPDF', 'wgtTMDPDF')
    h            : hadron number
    includeGluon : also fill the gluon (flavor 0) entry
    Xrange       : if omitted, defaults to X_RANGE_FF_DEFAULT for PDF='uTMDFF',
                   X_RANGE_PDF_DEFAULT otherwise
    """
    t0     = time.time()
    Xrange = _default_Xrange(Xrange, PDF)
    getter = _getter(PDF, kT=False)

    def TMDval_of(Qval, xval, bval):
        if xval > 0.999:
            return [0.] * 11
        return [xval * v for v in getter(xval, bval, h, Qval, Qval ** 2, includeGluon=includeGluon)]

    fill_flavors = _FLAVORS if includeGluon else _FLAVORS_NOG
    valuesList   = _build_grid([Qrange, Xrange, Brange], TMDval_of, _FLAVORS, fill_flavors)
    _write(path, [("Qg", Qrange), ("xg", Xrange), ("bg", Brange)], valuesList, t0)


def save_grid_kT(path, Qrange=Q_RANGE_DEFAULT, Xrange=None, Rrange=R_RANGE_DEFAULT,
                  PDF="uTMDPDF", h=1):
    """
    Save the TMDPDF/TMDFF grid in kT-space over a range of Q, Q^2, as returned
    by harpy (artemide must already be set up).

    PDF    : which TMD to save ('uTMDPDF', 'uTMDFF', 'SiversTMDPDF', 'wgtTMDPDF')
    h      : hadron number
    Xrange : if omitted, defaults to X_RANGE_FF_DEFAULT for PDF='uTMDFF',
             X_RANGE_PDF_DEFAULT otherwise

    Rrange is R = qT/Q; the physical kT = R*Q is passed to harpy.
    """
    t0     = time.time()
    Xrange = _default_Xrange(Xrange, PDF)
    getter = _getter(PDF, kT=True)

    def TMDval_of(Qval, xval, rval):
        if xval > 0.999:
            return [0.] * 11
        return [xval * v for v in getter(xval, rval * Qval, h, Qval, Qval ** 2, includeGluon=False)]

    valuesList = _build_grid([Qrange, Xrange, Rrange], TMDval_of, _FLAVORS_NOG, _FLAVORS_NOG)
    _write(path, [("Qg", Qrange), ("xg", Xrange), ("qToQg", Rrange)], valuesList, t0)
