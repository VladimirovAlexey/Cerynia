"""
Parsing of Tevatron W -> e nu transverse-momentum spectrum (D0, CDF).

Source: /data/arTeMiDe_Repository/data/CDF_D0/D0_W-98.json, CDF_W-91.json  (HEPData JSON)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY_W/*.csv

Ported from DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(W-Tevatron_old).py
-- definitions (process code, s, Q-range, y-range, thFactor, cuts, error
split) are kept identical to the old parsing.

NOTE: proc_id=6 here (W production), distinct from proc_id=3 (Z/gamma* near
peak) and proc_id=1 (generic/continuum virtual-photon DY) used elsewhere --
see project_nuclear_target_codes memory for the 1-vs-3 convention; 6 is a
new value for this category, not yet covered by that note.

Datasets renamed from the old "D0run1-W"/"CDFrun1-W" (and this script's own
"TevatronW_run1.py" filename) to drop the "run1" label and match the DY
branch's naming for the same Tevatron run 1 experiments (D01, CDF1), per
user instruction -- these are the same experiments, just the W-boson (not
Z/gamma*) final state:
    D01_W  -- D0 1.8 TeV, arXiv:hep-ph/9803003. Q-window is qT-dependent
              (kinematic bound M^2 > E1+E2 from the old code's own comment);
              rapidity restricted via y_min/y_max directly (not cutParams)
              since the source only restricts electron rapidity, per old
              parsing's comment.
    CDF1_W -- CDF 1.8 TeV, Phys.Rev.Lett. 66 (1991). qT bins are not given
              directly -- restored via a hardcoded variable-width binning
              function around each point's quoted centroid, ported verbatim
              as _cdf_qt_bin(); not re-derivable from the data alone
              ("high-qT bins assumed from the plot", per old parsing's own
              comment). y is fully unbounded here (rapidity restriction
              instead placed in cutParams), unlike D0 above -- inconsistent
              between the two datasets, but that is how old parsing had it;
              ported as-is.
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

import json
from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY_W/"

M_W = 80.  # forces <Q> = W-boson mass, matching the old parsing's explicit p["<Q>"]=M_W override

#%%
# ============================================================================
# D0 Run 1 -- arXiv:hep-ph/9803003
# ============================================================================
with open(path_to_data + "CDF_D0/D0_W-98.json") as f:
    values = json.load(f)["values"]

### Process is W-DY for p-pbar
ps_def, h_1, h_2, proc_id = 1, 1, -1, 6
s = 1800.**2
### rapidity restriction placed directly on y (not cutParams): "since they put
### restriction only on electron rapidity I do not restrict rapidity at all
### (it is done by y)" -- as in old parsing
y_min, y_max = -1.1, 1.1
includeCuts = True
cutParams = [25., 25., -20.1, 20.1]
### no luminosity/normalization error quoted (old parsing carried the append() commented out)
normErr = []

ds = DataSet.empty("DY", name="D01_W", comment="D0 1.8TeV normalized, W-boson production to electron",
                    reference="arXiv:hep-ph/9803003", normErr=normErr, isNormalized=True)

for i, v in enumerate(values):
    qT_min, qT_max = float(v["x"][0]["low"]), float(v["x"][0]["high"])
    y0 = v["y"][0]
    xSec = float(y0["value"])
    errs = [float(e["symerror"]) for e in y0["errors"]]  # stat, sys_1, sys_2

    ### Q-window is qT-dependent: lower edge from the kinematic bound M^2>E1+E2
    ### (old comment), floored at 10 GeV; upper edge fixed at 300 GeV (checked
    ### sufficient, per old parsing's comment). As in old parsing.
    Q_min = sqrt(50.**2 - qT_max**2) if qT_max < 50. else 2.
    if Q_min < 10.:
        Q_min = 10.
    Q_max = 300.

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"D01_W.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_W, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=errs[0],
        uncorrErr_1=errs[1],
        uncorrErr_2=errs[2],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")


#%%
# ============================================================================
# CDF Run 1 -- Phys.Rev.Lett. 66 (1991)
# ============================================================================
def _cdf_qt_bin(x):
    """Hardcoded variable-width qT binning around a point centroid (as in old parsing)."""
    if x < 20:
        return [x - 1., x + 1.]
    elif x < 50:
        return [x - 2.5, x + 2.5]
    elif x < 60:  # high-qT bins assumed from the plot (old parsing's own comment)
        return [50., 60.]
    elif x < 80:
        return [60., 80.]
    elif x < 130:
        return [80., 130.]
    else:
        return [130., 200.]


with open(path_to_data + "CDF_D0/CDF_W-91.json") as f:
    values = json.load(f)["values"]

ps_def, h_1, h_2, proc_id = 1, 1, -1, 6
s = 1800.**2
### no restriction on y (rapidity restriction instead placed in cutParams below), as in old parsing
y_min, y_max = -100., 100.
includeCuts = False
cutParams = [20., 20., -1.1, 1.1]
### fixed Q-window, lower boundary not defined in the source (as in old parsing)
Q_min, Q_max = 40., 300.
### no luminosity/normalization error quoted (old parsing carried the append() commented out)
normErr = []

ds = DataSet.empty("DY", name="CDF1_W", comment="CDF 1.8TeV, W-boson production to electron",
                    reference="Phys.Rev.Lett. 66 (1991)", normErr=normErr, isNormalized=True)

for i, v in enumerate(values):
    qT_centroid = float(v["x"][0]["value"])
    qT_min, qT_max = _cdf_qt_bin(qT_centroid)
    y0 = v["y"][0]
    xSec = float(y0["value"])
    errs = [float(e["symerror"]) for e in y0["errors"]]  # stat, sys

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CDF1_W.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_W, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max, qT_avg=qT_centroid,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=errs[0],
        uncorrErr_1=errs[1],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
