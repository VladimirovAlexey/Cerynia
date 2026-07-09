"""
Parsing of E537 (fixed-target, pbar+W -> mu+mu- Drell-Yan) invariant
cross-section data, Q-differential and xF-differential views.

Source: /data/arTeMiDe_Repository/data/FNAL-537/pbar+W(dQ).dat, pbar+W(dxF).dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "E537" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q/y-range, thFactor, error split) are kept identical to the
old parsing, including its hand-typed-table state-machine reader.

NOTE: thFactor is the fixed-target "invariant cross section" Jacobian,
1/(qT_max**2-qT_min**2)/(range_width)*0.001 -- structurally different from
the collider bin-integrated k/(qT_max-qT_min) form. Per user confirmation,
atmdeFactor still applies here: the old expression is multiplied by
atmdeFactor=(Q_max-Q_min)*(y_max-y_min)*(qT_max-qT_min), same as every other
DY case. Checked against DataProcessor.git's ReadDYdataFiles(low-energy).py:
thFactor is identical there for both dQ and dxF, so no ffactor-style fix
applies here (unlike E772/E605).

CORRECTED: dxF's h_1/h_2 order was swapped vs dQ in the originally-ported
main-file parsing ([2,1,-1,103] vs dQ's [2,-1,1,103]). The standalone
ReadDYdataFiles(low-energy).py uses dQ's order for both, with no swap --
taken as authoritative and applied here (h_1=-1 pbar / h_2=W target in both
datasets now). See project_e_series_thfactor_fix memory for the audit trail.

Raw format: comma-separated, hand-typed from the paper (not a HEPData export).
Q or xF sub-range markers appear as inline comment lines
("#: M(P=3 4) [GEV],,,4.0TO4.5" / "#: XL(...),,,-0.1TO0."); other "#:" lines,
blank lines, and the quoted column-header line are skipped. The raw table
gives PT**2 (center, low, high), not PT, so qT bounds are sqrt() of those.
Rows with zero-or-negative symmetrized uncorrErr (placeholder "not measured"
rows) are dropped, exactly as old parsing did.

Datasets (Phys.Rev.D 93 (1988) 1377):
    E537-dQ  -- Q-differential view (fixed y-window, Q swept per block)
    E537-dxF -- xF-differential view (fixed Q-window, y swept per block)
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from math import sqrt
from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"


def _read_e537_table(path, n_header, trigger_prefix):
    """
    Read a hand-typed FNAL-537 comma-separated table with inline sub-range
    marker lines ("#: <trigger_prefix>...<low>TO<high>"). Returns a list of
    (label, k, range_min, range_max, pt2_center, pt2_low, pt2_high, xSec, errp, errm)
    tuples, one per real data row -- label/range carried from the most recent
    marker line, matching the old state-machine parsing exactly. k is a
    per-label row index that resets to 0 at each new marker line and
    increments for every row seen there (whether or not it is later kept) --
    this matches the old parsing's point-id numbering exactly, including
    the "gaps" left by rows dropped for zero/negative uncorrErr.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]

    out = []
    label, rng, k = "?", (None, None), 0
    for line in lines[n_header:]:
        if line[:len(trigger_prefix)] == trigger_prefix:
            raw = line[-8:]
            label = raw.replace(".", "")
            lo, hi = raw.split("TO")
            rng = (float(lo), float(hi))
            k = 0
        elif line[:2] == "#:" or line == "" or line[:3] == '"PT':
            continue
        else:
            vals = [float(v) for v in line.split(",")]
            out.append((label, k, rng[0], rng[1], vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]))
            k += 1
    return out


#%%
# ============================================================================
# E537, Q-differential -- Phys.Rev.D 93 (1988) 1377
# ============================================================================
rows = _read_e537_table(path_to_data + "FNAL-537/pbar+W(dQ).dat", n_header=9, trigger_prefix="#: M")

### s_current as coded in old parsing (not a sqrt(s)**2 pattern like the other DY files -- ported verbatim)
s = 235.4
### Process is virtual-photon DY for pbar-tungsten
## ps_def, h_1, h_2, proc_id = 2, -1, 1, 103  ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 2, -1, 184074, 1  ##(>6.04 definition)
### fixed y-window for the Q-differential view (as in old parsing)
y_min, y_max = -0.1, 1.0
### 8% normalization uncertainty (as in old parsing)
normErr = [0.08]

ds = DataSet.empty("DY", name="E537-dQ", comment="E537 data Q-differential", reference="Phys.Rev.D 93 (1988) 1377",
                    normErr=normErr, isNormalized=False)

for label, k, Q_min, Q_max, pt2, pt2_min, pt2_max, xSec, errp, errm in rows:
    qT_min, qT_max = sqrt(pt2_min), sqrt(pt2_max)  ### table gives PT**2, not PT
    uncorrErr_0 = (errp - errm) / 2.
    if uncorrErr_0 <= 0:  # placeholder "not measured" rows (as in old parsing)
        continue

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"E537-dQ.{label}.{k}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### invariant cross-section Jacobian, divided by the (swept) Q window (as in old parsing)
        thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (Q_max - Q_min) * 0.001,
        includeCuts=False,
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# E537, xF-differential -- Phys.Rev.D 93 (1988) 1377
# ============================================================================
rows = _read_e537_table(path_to_data + "FNAL-537/pbar+W(dxF).dat", n_header=9, trigger_prefix="#: XL")

s = 235.4
### h_1/h_2 order fixed to match dQ (beam=pbar, target=W) -- the old main-file
### parsing had this swapped ([2,1,-1,103]) vs the standalone
### ReadDYdataFiles(low-energy).py, which uses the same order as dQ; the
### latter is treated as authoritative (see project_e_series_thfactor_fix)
##ps_def, h_1, h_2, proc_id = 2, -1, 1, 103  ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 2, -1, 184074, 1  ##(>6.04 definition)

### fixed Q-window for the xF-differential view (as in old parsing)
Q_min, Q_max = 4.0, 9.0
normErr = [0.08]

ds = DataSet.empty("DY", name="E537-dxF", comment="E537 data xF-differential", reference="Phys.Rev.D 93 (1988) 1377",
                    normErr=normErr, isNormalized=False)

for label, k, y_min, y_max, pt2, pt2_min, pt2_max, xSec, errp, errm in rows:
    qT_min, qT_max = sqrt(pt2_min), sqrt(pt2_max)
    uncorrErr_0 = (errp - errm) / 2.
    if uncorrErr_0 <= 0:
        continue

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"E537-dxF.{label}.{k}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### invariant cross-section Jacobian, divided by the (swept) y window (as in old parsing)
        thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * 0.001,
        includeCuts=False,
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
