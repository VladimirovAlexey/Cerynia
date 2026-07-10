"""
Parsing of JLab Hall A g2 (transverse spin structure function) data: a 2004
neutron measurement and two 2016 effective-neutron (He3 target) beam-energy
sub-tables.

Source: /data/arTeMiDe_Repository/data/g2Tables/JLab_Hall_A/2004/Table6.csv
        /data/arTeMiDe_Repository/data/g2Tables/JLab_Hall_A/2016/Tables_JLab_2016.csv
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/Hall_A_2004_neutron.csv,
        Hall_A(4.74 GeV).csv, Hall_A(5.89 GeV).csv (WW term)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "Hall A 2004"/"Hall A 2016" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins,
WW-term subtraction, and error treatment kept identical to the old parsing
(the user confirmed this script's physics/process definitions are already
in the modern standard; this is a pure reformatting port, nothing changed).

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction); proc_id
carried over unchanged from the old single `process` integer -- all THREE
datasets here use proc_id=101 (neutron), including the two "He3"-named
2016 sets: the old script's own inline comment on those lines says
"# it is neutron!" -- these are per-nucleon effective-neutron extractions
from a He3 target, not a distinct nuclear proc_id (matches the established
per-nucleon convention used elsewhere in this migration, e.g. deuteron
SIDIS targets).

NOTE: the old script also has a dead/fully-commented-out "Hall_A_04_He3"
block (2004 He3, not just neutron) -- never active, not ported.

NOTE: the 2004 dataset has real x_min/x_max AND real Q bin-guessing (both
guessed, via _guess_bin_sqrt/_guess_bin_linear, same Bottom/Middle/Top
logic as G2/E142.py/E155.py). The two 2016 sets are similarly both-guessed.

NOTE: the 2016 sets' raw table gives single-valued stat/syst errors (not
+/- pairs) -- used directly, NOT symmetrized via (a-b)/2 like every other
G2 dataset in this category (confirmed: old code uses `data[i][3]`,
`data[i][4]` raw, matches the old CSV's "Number of uncorr.errors,3" with
unsymmetrized values).

NOTE: normErr=[] throughout (matches every other G2 dataset).

Datasets:
    HallA-2004.n     -- reference nucl-ex/0405006v5, 3 points
    HallA-2016-4.He3 -- reference 1603.03612v3, E=4.74 GeV, 13 points, proc_id=101
    HallA-2016-5.He3 -- reference 1603.03612v3, E=5.89 GeV, 13 points, proc_id=101
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


def _guess_bin_linear(data, i, col):
    """Ported verbatim from old parsing's GuessBin_Bottom/Middle/Top("x",...) -- same logic, no sqrt."""
    val = lambda k: data[k][col]
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


#%% -- HallA-2004.n
### columns: [0]=x_avg,[1]=Q2_avg,[2]=g2(raw),[3]=stat+,[4]=stat-,[5]=syst+,[6]=syst-
kin = _read_csv_rows(path_to_data + "JLab_Hall_A/2004/Table6.csv", 21, 3)
ww = _read_csv_rows(path_to_ww + "Hall_A_2004_neutron.csv", 0, 3)
Qbins = [_guess_bin_sqrt(kin, i, 1) for i in range(len(kin))]
xbins = [_guess_bin_linear(kin, i, 0) for i in range(len(kin))]
s = 2 * 5.7 * M_neutron + M_neutron**2

ds = DataSet.empty("G2", name="HallA-2004.n", comment="Taken from table 6. Bin in Q guessed by us.",
                    reference="nucl-ex/0405006v5", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=101,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[1]),
        x_min=xbins[i][0], x_max=xbins[i][1], x_avg=row[0],
        thFactor=1.,
        xSec=row[2] - ww[i][0],
        uncorrErr_0=(row[3] - row[4]) / 2.,
        uncorrErr_1=(row[5] - row[6]) / 2.,
        uncorrErr_2=ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- HallA-2016-4.He3 (E=4.74 GeV; process is neutron, "it is neutron!" per old parsing's own comment)
### columns: [0]=x_avg,[1]=Q2_avg,[2]=g2(raw),[3]=stat(raw, unsymmetrized),[4]=syst(raw, unsymmetrized)
kin = _read_csv_rows(path_to_data + "JLab_Hall_A/2016/Tables_JLab_2016.csv", 2, 13)
ww = _read_csv_rows(path_to_ww + "Hall_A(4.74 GeV).csv", 0, 13)
Qbins = [_guess_bin_sqrt(kin, i, 1) for i in range(len(kin))]
xbins = [_guess_bin_linear(kin, i, 0) for i in range(len(kin))]
s = 2 * 4.74 * M_neutron + M_neutron**2

ds = DataSet.empty("G2", name="HallA-2016-4.He3", comment="E=4.74 GeV. Taken from table VIII in paper. Bin in Q guessed by us.",
                    reference="1603.03612v3", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=101,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[1]),
        x_min=xbins[i][0], x_max=xbins[i][1], x_avg=row[0],
        thFactor=1.,
        xSec=row[2] - ww[i][0],
        uncorrErr_0=row[3],
        uncorrErr_1=row[4],
        uncorrErr_2=ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- HallA-2016-5.He3 (E=5.89 GeV; process is neutron, "it is neutron!" per old parsing's own comment)
kin = _read_csv_rows(path_to_data + "JLab_Hall_A/2016/Tables_JLab_2016.csv", 16, 13)
ww = _read_csv_rows(path_to_ww + "Hall_A(5.89 GeV).csv", 0, 13)
Qbins = [_guess_bin_sqrt(kin, i, 1) for i in range(len(kin))]
xbins = [_guess_bin_linear(kin, i, 0) for i in range(len(kin))]
s = 2 * 5.89 * M_neutron + M_neutron**2

ds = DataSet.empty("G2", name="HallA-2016-5.He3", comment="E=5.89 GeV. Taken from table IX in paper. Bin in Q guessed by us.",
                    reference="1603.03612v3", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=101,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[1]),
        x_min=xbins[i][0], x_max=xbins[i][1], x_avg=row[0],
        thFactor=1.,
        xSec=row[2] - ww[i][0],
        uncorrErr_0=row[3],
        uncorrErr_1=row[4],
        uncorrErr_2=ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")
