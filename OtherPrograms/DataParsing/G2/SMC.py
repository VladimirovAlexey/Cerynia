"""
Parsing of SMC g2 (transverse spin structure function) data, proton.

Source: /data/arTeMiDe_Repository/data/g2Tables/SMC/Table8.csv, Table9.csv
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "SMC" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins, and
error treatment kept identical to the old parsing (the user confirmed this
script's physics/process definitions are already in the modern standard;
this is a pure reformatting port, nothing changed).

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction); proc_id
carried over unchanged from the old single `process` integer (100=proton).

NOTE: unlike every other G2 dataset in this category, SMC does NOT need a
separately-computed WW-term subtraction: Table8 gives the raw g2, Table9
gives "G2 (C=WW)" (already WW-corrected by the experiment itself, per the
raw file's own column header) -- xSec is the direct difference between the
two tables' values at the same row index (assumed row-for-row
correspondence between Table8 and Table9, per old parsing). Systematic
uncertainty is neglected/set to 0 by the experiment (per old comment), so
only 2 uncorrErr columns (Table8's symmetrized stat error, Table9's
symmetrized stat error).

NOTE: only Table8's kinematics (x/Q bins, both guessed via
_guess_bin_sqrt/_guess_bin_linear, same Bottom/Middle/Top logic as
G2/E142.py) are used for the point's x/Q_min/Q_max/Q_avg/x_avg; Table9
contributes only its own value/error columns.

NOTE: normErr=[] (matches every other G2 dataset).

Dataset (arXiv:hep-ex/9702005), 6 points:
    SMC.p
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/g2Tables/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/G2/"

M_proton = 0.938

ps_def, h_1 = 1, 1


def _guess_bin_sqrt(data, i, col):
    """Ported verbatim from old parsing's GuessBin_Bottom/Middle/Top("Q",...), see G2/E142.py."""
    val = lambda k: sqrt(data[k][col])
    if i == 0:
        return [val(i) - (val(i + 1) - val(i)) / 2, val(i) + (val(i + 1) - val(i)) / 2]
    elif i == len(data) - 1:
        return [val(i) - (val(i) - val(i - 1)) / 2, val(i) + (val(i) - val(i - 1)) / 2]
    else:
        return [val(i) - (val(i) - val(i - 1)) / 2, val(i) + (val(i + 1) - val(i)) / 2]


def _read_csv_rows(path, start, n_rows):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    return [[float(v) for v in line.split(",")] for line in lines[start:start + n_rows]]


#%% -- proton
### columns: [0]=x_avg,[1]=x_min,[2]=x_max,[3]=Q2_avg,[4]=y(unused),[5]=value,[6]=stat+,[7]=stat-
kin1 = _read_csv_rows(path_to_data + "SMC/Table8.csv", 13, 6)
kin2 = _read_csv_rows(path_to_data + "SMC/Table9.csv", 13, 6)
Qbins = [_guess_bin_sqrt(kin1, i, 3) for i in range(len(kin1))]

s = 2 * 190 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="SMC.p",
                    comment="Taken from tables 8 and 9. Bin in Q guessed by us. They neglect the systematic uncertainty on g2 and setting it to 0. Furthermore, they compute the WW-terms and add their uncertainty.",
                    reference="arXiv:hep-ex/9702005", normErr=[], isNormalized=False)

for i in range(len(kin1)):
    row1, row2 = kin1[i], kin2[i]
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=100,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row1[3]),
        x_min=row1[1], x_max=row1[2], x_avg=row1[0],
        thFactor=1.,
        xSec=row1[5] - row2[5],
        uncorrErr_0=(row1[6] - row1[7]) / 2.,
        uncorrErr_1=(row2[6] - row2[7]) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
