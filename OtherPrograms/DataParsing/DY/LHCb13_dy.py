"""
Parsing of LHCb 13 TeV Z/gamma* -> mu+mu- transverse-momentum spectrum,
2021 update, double-differential in (y, qT).

Source: /data/arTeMiDe_Repository/data/LHCb/LHCb_13_dy[2021].dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/LHCb13_dy.csv

Ported from the "LHCb13 y-differential" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(LHCb13).py --
definitions (process code, s, Q-range, thFactor, cuts, error split) are kept
identical to the old parsing. Only the "LHCb13_dy" case from that file is
handled here (the y-integrated "LHCb13(2021)" case is a separate, not yet
ported, dataset). Dataset renamed from the old "LHCb13_dy(2021)" to
"LHCb13_dy" per user instruction (drop the "(2021)" label).

Raw format: hand-typed from arXiv:2112.07458, comma-separated, 6 header lines
+ 1 comment line (7 total), 70 data rows (5 y-bins x 14 qT-bins), no trailing
blank line. Columns [yMin, yMax, ptMin, ptMax, xSec, stat, sys, lumi(2%)] --
the lumi column is redundant with the global 2% normErr and unused per-point,
exactly as old parsing. The systematic (7th column) is lightly correlated and
kept uncorrelated (ignored as fully correlated), exactly as old parsing.

Dataset:
    LHCb13_dy -- arXiv:2112.07458
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.1876  # this script's own value (matches CMS13.py), ported as coded

#%%
# ============================================================================
# LHCb 13 TeV, y-differential, 2021 update -- arXiv:2112.07458
# ============================================================================
with open(path_to_data + "LHCb/LHCb_13_dy[2021].dat") as f:
    lines = [line.rstrip("\n") for line in f]
rows = [[float(v) for v in line.split(",")] for line in lines[7:]]

### M(mu mu) window and PT(mu)>20 GeV cut, stated in the file header
Q_min, Q_max = 60., 120.
### 13 TeV
s = 13000.**2
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### ETARAP(MU) : 2.0 - 4.5, stated in the file header (as fiducial cut)
includeCuts = True
cutParams = [20., 20., 2., 4.5]
### 2% luminosity uncertainty (table 1)
normErr = [0.02]

ds = DataSet.empty("DY", name="LHCb13_dy", comment="LHCb 13TeV 2021 update y-differential",
                    reference="arXiv:2112.07458", normErr=normErr, isNormalized=False)

for i, (y_min, y_max, qT_min, qT_max, xSec, stat, sys, lumi) in enumerate(rows):
    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"LHCb13_dy.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size in both qT and y (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min) / (y_max - y_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        ### stat + syst, both treated uncorrelated (syst lightly correlated, ignored as in old parsing)
        uncorrErr_0=stat,
        uncorrErr_1=sys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
