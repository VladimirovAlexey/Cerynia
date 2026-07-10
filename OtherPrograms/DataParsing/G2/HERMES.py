"""
Parsing of HERMES g2 (transverse spin structure function) data, proton:
per-(x,Q^2)-bin table and the Q^2-evolved/x-averaged companion table.

Source: /data/arTeMiDe_Repository/data/g2Tables/HERMES/Table1.csv, Table2.csv
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/HERMES.csv, HERMES_AVG.csv (WW term)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "HERMES"/"HERMES.av" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins,
WW-term subtraction, and error treatment kept identical to the old parsing
(the user confirmed this script's physics/process definitions are already
in the modern standard; this is a pure reformatting port, nothing changed).

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction); proc_id
carried over unchanged from the old single `process` integer (100=proton
for both datasets).

NOTE: the raw table's "xg2" column is x*g2, not g2 -- xSec = raw_xg2 -
x*WW_value, thFactor=x (=<x>), same "x-scaled" pattern as some E155/E143
sets. x bins come directly from the raw table's own columns (real, not
guessed).

NOTE: the Q window is a FIXED literal [sqrt(0.18), sqrt(20.)] for every
point in both datasets -- the old script also computes a real per-point Q
bin guess (Qbounds1/2, via the same Bottom/Middle/Top logic used
elsewhere) but never actually uses it for `p["Q"]`; only `<Q>` is set from
the measured average. Confirmed dead code by cross-checking the old CSV
(every point has the same Qmin/Qmax, only <Q> varies) -- ported verbatim,
including not computing the unused Qbounds guess.

NOTE: normErr=[] (matches every other G2 dataset). The old dataset comment
references an external correlation-matrix file
("Correlation_Matrix_23BINS.csv"/"...9BINS.csv") that is never actually
loaded or used anywhere in the old script -- informational only, comment
ported verbatim.

Datasets (arXiv:1112.5584v3):
    HERMES     -- Table 1, per-(x,Q^2)-bin, 23 points
    HERMES.av  -- Table 2, Q^2-evolved/x-averaged, 9 points
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/g2Tables/"
path_to_ww = "/data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/G2/"

M_proton = 0.938

ps_def, h_1 = 1, 1
### fixed Q window (as in old parsing; the per-point Qbounds guess is computed but unused, see module NOTE)
Q_min, Q_max = sqrt(0.18), sqrt(20.)
s = 2 * 27.6 * M_proton + M_proton**2


def _read_csv_rows(path, start, n_rows):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    return [[float(v) for v in line.split(",")] for line in lines[start:start + n_rows]]


### columns: [0]=BIN(unused),[1]=X-RANGE(unused),[2]=x_min,[3]=x_max,[4]=x_avg,[5]=Q2_avg,[6]=xg2(raw),[7]=stat+,[8]=stat-,[9]=syst+,[10]=syst-

#%% -- HERMES (Table 1, per-(x,Q^2) bin)
kin = _read_csv_rows(path_to_data + "HERMES/Table1.csv", 12, 23)
ww = _read_csv_rows(path_to_ww + "HERMES.csv", 0, 23)

ds = DataSet.empty("G2", name="HERMES",
                    comment="Taken from table 1. Bin in Q guessed by us. Errors in g2 are correlated between bins over x. The correlation matrix is saved in: '/home/guillermo/Work/Twist_3/Exp_Data/HERMES/Correlation_Matrix_23BINS.csv' ",
                    reference="arXiv:1112.5584v3", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    x_avg = row[4]
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=100,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[5]),
        x_min=row[2], x_max=row[3], x_avg=x_avg,
        thFactor=x_avg,
        xSec=row[6] - x_avg * ww[i][0],
        uncorrErr_0=(row[7] - row[8]) / 2.,
        uncorrErr_1=(row[9] - row[10]) / 2.,
        uncorrErr_2=x_avg * ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- HERMES.av (Table 2, Q^2-evolved/x-averaged)
kin = _read_csv_rows(path_to_data + "HERMES/Table2.csv", 12, 9)
ww = _read_csv_rows(path_to_ww + "HERMES_AVG.csv", 0, 9)

ds = DataSet.empty("G2", name="HERMES.av",
                    comment="Taken from table 2. Bin in Q guessed by us. Errors in g2 are correlated between bins over x. The correlation matrix is saved in: '/home/guillermo/Work/Twist_3/Exp_Data/HERMES/Correlation_Matrix_9BINS.csv' ",
                    reference="arXiv:1112.5584v3", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    x_avg = row[4]
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=100,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[5]),
        x_min=row[2], x_max=row[3], x_avg=x_avg,
        thFactor=x_avg,
        xSec=row[6] - x_avg * ww[i][0],
        uncorrErr_0=(row[7] - row[8]) / 2.,
        uncorrErr_1=(row[9] - row[10]) / 2.,
        uncorrErr_2=x_avg * ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")
