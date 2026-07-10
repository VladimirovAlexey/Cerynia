"""
Parsing of E143 g2 (transverse spin structure function) data: two separate
E143 sub-measurements -- the 1998 g2-bar (higher-twist) paper (p/d/n) and
the 1996 raw g2 HEPData release (p/d, needing WW-term subtraction).

Source: /data/arTeMiDe_Repository/data/g2Tables/E143/1998/*.csv (g2-bar)
        /data/arTeMiDe_Repository/data/g2Tables/E143/1996/*.csv (raw g2)
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/E143_1996_*.csv (WW term, 1996 sub-case only)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "E143 1998 proton/deuteron/neutron" and "E143 1995/1996
proton/deuteron" blocks of DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py
-- s, Q/x bins, WW-term subtraction, and error treatment kept identical to
the old parsing (the user confirmed this script's physics/process
definitions are already in the modern standard; this is a pure
reformatting port, nothing changed).

NOTE: process code. Old parsing used a single `process` integer (100=p,
101=n, 102=d). Per user instruction, this maps directly onto Cerynia's G2
schema as proc_id (unchanged value), with ps_def=1, h_1=1 fixed.

NOTE: the 1998 g2-bar sets (E143.p/d/n) are the ONLY G2 sub-case in this
category presenting g2-bar directly (higher-twist part with WW already
removed by the experiment itself) -- xSec is read straight from the
twist3/g2bar data file, no WW subtraction needed here (unlike every other
G2 dataset). Only 2 uncorrErr columns (stat, syst of g2-bar itself).

NOTE: the 1995/1996 sets (E143-1995.p/d) DO need WW subtraction:
xSec = raw_g2 - WW_value, uncorrErr = [stat(symmetrized),
syst(symmetrized), WW_stat] (3 columns). Each of these two datasets merges
TWO independently-Q-bin-guessed raw tables (Table1+Table2 for proton,
Table3+Table4 for deuteron) end to end.

NOTE: E143-1995.d uses M_deuteron for s (unlike E143.d in the 1998 set,
which uses M_proton for s) -- ported verbatim, not harmonized.

NOTE: Q bin is "guessed" via the old script's GuessBin_Bottom/Middle/Top
helpers (see _guess_bin_sqrt below, same as G2/E142.py). For E143.p/d/n
the 12 rows form TWO independent runs (rows 0-6 and 7-11, a discontinuity
in the raw table between two spectrometer settings) -- each run gets its
own Bottom/Middle/Top guess, not one continuous 12-row run. For
E143-1995.p/d each of the two merged tables (7 rows + 5 rows) is likewise
its own independent run.

NOTE: normErr=[] throughout -- the old script computes `lumUncertainty`
per dataset (3.7%/4.9%/4.9% for the 1998 p/d/n sets) but never appends it
to any DataSet's normErr (confirmed by every old G2 CSV: "Number of
norm.errors,0"). Also note the old comments for E143-1995.p/d say "Taken
from tables 1,2, 3" but the code only reads 2 tables each (a leftover from
copy-pasting) -- ported verbatim, not corrected.

Datasets:
    E143.p, E143.d, E143.n       -- reference hep-ph/9802357, 12 points each
    E143-1995.p, E143-1995.d     -- reference 10.17182/hepdata.19584.v1, 12 points each
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/g2Tables/"
path_to_ww = "/data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/G2/"

M_proton = 0.938
M_neutron = 0.939
M_deuteron = (M_proton + M_neutron) / 2

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


def _qbin_two_runs(run1, run2, col):
    """Independent Bottom/Middle/Top guess per run, then concatenated (as in old parsing's split-table datasets)."""
    return [_guess_bin_sqrt(run1, i, col) for i in range(len(run1))] + \
           [_guess_bin_sqrt(run2, i, col) for i in range(len(run2))]


#%% -- E143.p (1998 g2-bar)
### columns: [0]=x_avg,[1]=x_min,[2]=x_max,[3]=Q2_avg (only cols 0-3 used)
kin = _read_csv_rows(path_to_data + "E143/1998/Table31.csv", 14, 12)
### twist3/g2bar columns: [0]=xSec(g2bar), [1]=stat, [2]=syst
g2bar = _read_csv_rows(path_to_data + "E143/1998/E143_1998_g2bar_proton.csv", 0, 12)
Qbins = _qbin_two_runs(kin[0:7], kin[7:12], 3)
s = 2 * 29.1 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="E143.p",
                    comment="g2bar data taken from page 86 in the paper. Bin in Q guessed by us. Additional normalization uncertainty 3.7%. The E143 experiment is the only one which presents g2-bar, therefore the uncertainties are solely the stat+syst of said g2-bar, unlike the rest of the experimental tables",
                    reference="hep-ph/9802357", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=100,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=g2bar[i][0],
        uncorrErr_0=g2bar[i][1],
        uncorrErr_1=g2bar[i][2],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E143.d (1998 g2-bar)
kin = _read_csv_rows(path_to_data + "E143/1998/Table32.csv", 12, 12)
g2bar = _read_csv_rows(path_to_data + "E143/1998/E143_1998_g2bar_deuteron.csv", 0, 12)
Qbins = _qbin_two_runs(kin[0:7], kin[7:12], 3)
s = 2 * 29.1 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="E143.d",
                    comment="g2bar data taken from page 86 in the paper. Bin in Q guessed by us. Additional normalization uncertainty 4.9%. The E143 experiment is the only one which presents g2-bar, therefore the uncertainties are solely the stat+syst of said g2-bar, unlike the rest of the experimental tables",
                    reference="hep-ph/9802357", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=102,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=g2bar[i][0],
        uncorrErr_0=g2bar[i][1],
        uncorrErr_1=g2bar[i][2],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E143.n (1998 g2-bar)
kin = _read_csv_rows(path_to_data + "E143/1998/Table33.csv", 14, 12)
g2bar = _read_csv_rows(path_to_data + "E143/1998/E143_1998_g2bar_neutron.csv", 0, 12)
Qbins = _qbin_two_runs(kin[0:7], kin[7:12], 3)
s = 2 * 29.1 * M_neutron + M_neutron**2

ds = DataSet.empty("G2", name="E143.n",
                    comment="g2bar data taken from page 87 in the paper. Bin in Q guessed by us. Additional normalization uncertainty 4.9%. The E143 experiment is the only one which presents g2-bar, therefore the uncertainties are solely the stat+syst of said g2-bar, unlike the rest of the experimental tables",
                    reference="hep-ph/9802357", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=101,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=g2bar[i][0],
        uncorrErr_0=g2bar[i][1],
        uncorrErr_1=g2bar[i][2],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E143-1995.p (1996 raw g2, WW-subtracted)
### columns: [0]=x_avg,[1]=x_min,[2]=x_max,[3]=Q2_avg,[4]=g2(raw),[5]=stat+,[6]=stat-,[7]=syst+,[8]=syst-
kin1 = _read_csv_rows(path_to_data + "E143/1996/Table1.csv", 25, 7)
kin2 = _read_csv_rows(path_to_data + "E143/1996/Table2.csv", 23, 5)
kin = kin1 + kin2
ww1 = _read_csv_rows(path_to_ww + "E143_1996_1P.csv", 0, 7)
ww2 = _read_csv_rows(path_to_ww + "E143_1996_2P.csv", 0, 5)
ww = ww1 + ww2
Qbins = _qbin_two_runs(kin1, kin2, 3)
s = 2 * 29.1 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="E143-1995.p",
                    comment="Taken from tables 1,2, 3. Bin in Q guessed by us.",
                    reference="10.17182/hepdata.19584.v1", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=100,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=row[4] - ww[i][0],
        uncorrErr_0=(row[5] - row[6]) / 2.,
        uncorrErr_1=(row[7] - row[8]) / 2.,
        uncorrErr_2=ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E143-1995.d (1996 raw g2, WW-subtracted)
kin1 = _read_csv_rows(path_to_data + "E143/1996/Table3.csv", 22, 7)
kin2 = _read_csv_rows(path_to_data + "E143/1996/Table4.csv", 20, 5)
kin = kin1 + kin2
ww1 = _read_csv_rows(path_to_ww + "E143_1996_1D.csv", 0, 7)
ww2 = _read_csv_rows(path_to_ww + "E143_1996_2D.csv", 0, 5)
ww = ww1 + ww2
Qbins = _qbin_two_runs(kin1, kin2, 3)
s = 2 * 29.1 * M_deuteron + M_deuteron**2

ds = DataSet.empty("G2", name="E143-1995.d",
                    comment="Taken from tables 1,2, 3. Bin in Q guessed by us.",
                    reference="10.17182/hepdata.19584.v1", normErr=[], isNormalized=False)
for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=102,
        s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=row[4] - ww[i][0],
        uncorrErr_0=(row[5] - row[6]) / 2.,
        uncorrErr_1=(row[7] - row[8]) / 2.,
        uncorrErr_2=ww[i][1],
    ))
ds.save_csv(path_to_save + ds.name + ".csv")
