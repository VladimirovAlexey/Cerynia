"""
Parsing of ATLAS 13 TeV Z/gamma* -> l+l- transverse-momentum spectrum,
normalized to 1/sigma, rapidity-inclusive.

Source: /data/arTeMiDe_Repository/data/ATLAS/A13.json  (HEPData JSON export)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/A13-norm.csv

Ported from the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(ATLAS13).py --
definitions (process code, s, Q/y-range, thFactor, cuts, error split) are
kept identical to the old parsing.

Raw format: HEPData JSON, "values" is a flat list of 44 points, each with
x[0].low/high (qT bin) and two y-entries: y[0] (group 0) is the actual
(1/sigma)dsigma/dqT value with two percentage errors (errors[0]=correlated,
errors[1]=uncorrelated, both "N.NNNN%" strings); y[1] (group 1) is an
unrelated correction/acceptance factor, unused, exactly as old parsing.

Dataset:
    A13-norm -- arXiv:1912.02844
"""

import sys
import json
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.  # forces <Q> = Z-boson mass rather than the (very close) bin midpoint,
           # matching the old parsing's explicit p["<Q>"]=M_Z override

#%%
# ============================================================================
# ATLAS 13 TeV -- arXiv:1912.02844
# ============================================================================
with open(path_to_data + "ATLAS/A13.json") as f:
    values = json.load(f)["values"]

### M(ll) window and lepton pT>27 GeV, |eta|<2.5 cuts (as in old parsing)
Q_min, Q_max = 66., 116.
### 13 TeV
s = 13000.**2
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### full rapidity window (as in old parsing; no symmetrize factor needed)
y_min, y_max = -2.5, 2.5
includeCuts = True
cutParams = [27., 27., -2.5, 2.5]
### normalized (1/sigma)dsigma/dqT measurement -- no luminosity/normalization
### error (old parsing's normErr.append(lumUncertainty) line was commented out)
normErr = []

ds = DataSet.empty("DY", name="A13-norm", comment="ATLAS 13TeV normalized to 1/sigma",
                    reference="arXiv:1912.02844", normErr=normErr, isNormalized=True)

for i, v in enumerate(values):
    qT_min, qT_max = float(v["x"][0]["low"]), float(v["x"][0]["high"])
    y0 = v["y"][0]  # group 0 = the actual cross-section value; y[1] (group 1) is unused, as in old parsing
    xSec = float(y0["value"])
    ### errors are given as percentages ("N.NNNN%"): errors[0]=correlated, errors[1]=uncorrelated
    corrErr_0 = xSec * 0.01 * float(y0["errors"][0]["symerror"][:-1])
    uncorrErr_0 = xSec * 0.01 * float(y0["errors"][1]["symerror"][:-1])

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A13-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size, no symmetrize factor (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=corrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
