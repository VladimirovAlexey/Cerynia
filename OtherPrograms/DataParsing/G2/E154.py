"""
Parsing of E154 g2 (transverse spin structure function) data, neutron.

Source: /data/arTeMiDe_Repository/data/g2Tables/E154/Table1.csv, Table2.csv
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/E154_1N.csv, E154_2N.csv (WW term)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "E154 neutron" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins,
WW-term subtraction, and error treatment kept identical to the old parsing
(the user confirmed this script's physics/process definitions are already
in the modern standard; this is a pure reformatting port, nothing changed).

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction); proc_id
carried over unchanged from the old single `process` integer (101=neutron).

NOTE: unlike most other G2 sub-cases, x bins here come directly from the
raw table's own x_min/x_max columns -- NOT guessed. Only the Q bin is
guessed (via _guess_bin_sqrt, same helper as G2/E142.py/E143.py), applied
independently to each of the two merged tables (10-row Table1 + 7-row
Table2 -- two separate spectrometer-setting runs, each with its own
Bottom/Middle/Top guess, not one continuous 17-row run).

NOTE: this dataset's old script block has a reference/comment that is
almost certainly a copy-paste leftover from the neighboring E155-1999
block ("hep-ex/9901006v1", "Taken from tables 1, 2, 3... No syst.
uncert.") -- it doesn't match E154's actual behavior (2 tables read, not
3; a syst error term IS included in uncorrErr, contradicting "No syst.
uncert."). Ported verbatim, not corrected -- flagged per the user's
"nothing to change" instruction.

NOTE: normErr=[] (matches every other G2 dataset; old script's
`lumUncertainty` is computed but never appended anywhere in this category).

Dataset, 17 points (10+7):
    E154.n
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/g2Tables/"
path_to_ww = "/data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/G2/"

M_neutron = 0.939

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


#%% -- neutron
### columns: [0]=x_avg,[1]=x_min,[2]=x_max,[3]=Q2_avg,[4]=g2(raw),[5]=stat+,[6]=stat-,[7]=syst+,[8]=syst-
kin1 = _read_csv_rows(path_to_data + "E154/Table1.csv", 28, 10)
kin2 = _read_csv_rows(path_to_data + "E154/Table2.csv", 25, 7)
kin = kin1 + kin2
ww1 = _read_csv_rows(path_to_ww + "E154_1N.csv", 0, 10)
ww2 = _read_csv_rows(path_to_ww + "E154_2N.csv", 0, 7)
ww = ww1 + ww2
Qbins = [_guess_bin_sqrt(kin1, i, 3) for i in range(len(kin1))] + \
        [_guess_bin_sqrt(kin2, i, 3) for i in range(len(kin2))]

s = 2 * 48.3 * M_neutron + M_neutron**2

ds = DataSet.empty("G2", name="E154.n",
                    comment="Taken from tables 1, 2, 3. Bin in Q and x guessed by us. No syst. uncert.",
                    reference="hep-ex/9901006v1", normErr=[], isNormalized=False)

for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=101,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=row[4] - ww[i][0],
        uncorrErr_0=(row[5] - row[6]) / 2.,
        uncorrErr_1=(row[7] - row[8]) / 2.,
        uncorrErr_2=ww[i][1],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
