"""
Parsing of CMS 8 TeV W -> l nu transverse-momentum spectrum (electron, muon).

Source: /data/arTeMiDe_Repository/data/CMS/WbosonELECTRON_2016.json, WbosonMUON_2016.json  (HEPData JSON)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY_W/*.csv

Ported from DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(W_CMS).py
-- definitions (process code, s, Q-range, y-range, thFactor, cuts, error
split) are kept identical to the old parsing.

NOTE: proc_id=6 (W production), same convention as TevatronW.py.
cutParams_1 (kCut2) is 0. for both datasets in old parsing -- W decays to a
single visible lepton + neutrino, so there is no second-lepton cut to encode
(unlike the DY branch's Z/gamma* cutParams, which cut on both leptons).

Datasets (arXiv:1606.05864):
    CMS8_W-electron -- 25 GeV lepton-pT cut, |y|<2.5
    CMS8_W-muon     -- 20 GeV lepton-pT cut, |y|<2.1 (note cutParams' eta
                       window is still -2.5/2.5 for this dataset too, not
                       -2.1/2.1 matching its own y -- ported verbatim from
                       old parsing, not corrected)
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

import json
from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY_W/"

M_W = 80.  # forces <Q> = W-boson mass, matching the old parsing's explicit p["<Q>"]=M_W override

#%%
# ============================================================================
# CMS 8 TeV, W -> e nu -- arXiv:1606.05864
# ============================================================================
with open(path_to_data + "CMS/WbosonELECTRON_2016.json") as f:
    values = json.load(f)["values"]

### Process is W-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 6
s = 8000.**2
Q_min, Q_max = 20., 300.
y_min, y_max = -2.5, 2.5
includeCuts = True
cutParams = [25., 0., -2.5, 2.5]
### no luminosity/normalization error quoted (old parsing carried the append() commented out)
normErr = []

ds = DataSet.empty("DY", name="CMS8_W-electron", comment="CMS 8TeV normalized, W-boson production to electron",
                    reference="arXiv:1606.05864", normErr=normErr, isNormalized=True)

for i, v in enumerate(values):
    qT_min, qT_max = float(v["x"][0]["low"]), float(v["x"][0]["high"])
    y0 = v["y"][0]
    xSec = float(y0["value"])
    uncorrErr_0 = float(y0["errors"][0]["symerror"])

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS8_W-electron.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_W, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# CMS 8 TeV, W -> mu nu -- arXiv:1606.05864
# ============================================================================
with open(path_to_data + "CMS/WbosonMUON_2016.json") as f:
    values = json.load(f)["values"]

ps_def, h_1, h_2, proc_id = 1, 1, 1, 6
s = 8000.**2
Q_min, Q_max = 20., 300.
y_min, y_max = -2.1, 2.1
includeCuts = True
### eta window here is -2.5/2.5, not this dataset's own -2.1/2.1 y-window -- ported verbatim from old parsing
cutParams = [20., 0., -2.5, 2.5]
normErr = []

ds = DataSet.empty("DY", name="CMS8_W-muon", comment="CMS 8TeV normalized, W-boson production to muon",
                    reference="arXiv:1606.05864", normErr=normErr, isNormalized=True)

for i, v in enumerate(values):
    qT_min, qT_max = float(v["x"][0]["low"]), float(v["x"][0]["high"])
    y0 = v["y"][0]
    xSec = float(y0["value"])
    uncorrErr_0 = float(y0["errors"][0]["symerror"])

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS8_W-muon.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_W, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
