"""
Parsing of E155 g2 (transverse spin structure function) data: three beam
energies (29.1, 32.3 GeV from the 2003 paper; 38.8 GeV from the 1999
paper), proton and deuteron each.

Source: /data/arTeMiDe_Repository/data/g2Tables/E155/2003/Table1-6.csv
        /data/arTeMiDe_Repository/data/g2Tables/E155/1999/Table1-3.csv
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/E155_*.csv (WW term)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "E155 2003"/"E155 1999" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins,
WW-term subtraction, systematic-error formulas, and error treatment kept
identical to the old parsing (the user confirmed this script's physics/
process definitions are already in the modern standard; this is a pure
reformatting port, nothing changed).

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction); proc_id
carried over unchanged from the old single `process` integer (100=proton,
102=deuteron).

NOTE: raw table columns here are [x, <Q^2>, x*g2(raw), error+, error-] --
note the raw "g2" column is actually x*g2 (confirmed by the raw CSV's own
column header "X*G2"), not g2 itself. thFactor is therefore set to x
(=data[i][0]) for most of these datasets, to convert a plain-g2 theory
prediction into the x*g2 form the data is given in -- EXCEPT it is
inconsistently 1 for some deuteron sets, ported exactly as coded per
dataset (see table below, verified against the old script line by line,
not assumed to be uniform):
    E155-29.p: thFactor=x   E155-29.d: thFactor=x
    E155-32.p: thFactor=x   E155-32.d: thFactor=1
    E155-38.p: thFactor=x   E155-38.d: thFactor=x

(verified directly against the old script's line numbers, not assumed
uniform -- an initial pass mis-extracted E155-29.d as thFactor=1; the old
script's actual line 1703 is `p["thFactor"]=data_current[i][0]`, caught
and fixed by cross-checking against the old CSV's own Th.Factor column)

NOTE: each dataset merges THREE raw tables end to end, each independently
Q-bin- and x-bin-guessed (via _guess_bin_sqrt/_guess_bin_linear, same
Bottom/Middle/Top logic as G2/E142.py, applied per-table not across the
merged whole).

NOTE: xSec = raw_x*g2 - x*WW_value (WW term scaled by x, unlike the
E142/E143-1995/E154 datasets which subtract WW unscaled). uncorrErr =
[error(symmetrized), syst(x) (from the formula printed in the raw file's
own "#: SYS," header line), x*WW_stat] for the 29.1/32.3 GeV sets (3
columns); the 38.8 GeV sets have NO syst term ("No syst. uncert." per the
raw file, no "#: SYS," line present) so only 2 uncorrErr columns
(error(symmetrized), x*WW_stat).

NOTE: E155-32 (32.3 GeV) reuses the SAME WW files as E155-29 (29.1 GeV) --
old script's own comment: "THE POINTS FOR x AND Q^2 ARE THE SAME FOR BOTH
29.1 and 32.3 GeV, THIS MEANS THAT WE ONLY NEED THE WW-TERMS FROM THE
PREVIOUS COMPUTATION". Ported verbatim (not re-derived).

NOTE: E155-32.d's old comment says "Taken from tables 1,2, 3" -- a
copy-paste leftover from the 29.1 GeV deuteron block (the 32.3 GeV block
actually reads Table4/5/6). Ported verbatim, not corrected.

NOTE: normErr=[] throughout (matches every other G2 dataset).

Datasets (arXiv:hep-ex/0204028v1 for 29.1/32.3 GeV, arXiv:hep-ex/9901006v1
for 38.8 GeV):
    E155-29.p, E155-29.d   -- 20 points each (8+7+5)
    E155-32.p, E155-32.d   -- 20 points each (8+7+5)
    E155-38.p, E155-38.d   -- 22 points each (10+7+5)
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


def _merge_bins(runs, col, guess_fn):
    out = []
    for run in runs:
        out += [guess_fn(run, i, col) for i in range(len(run))]
    return out


def syst_p(x):
    return 0.0016 - x * 0.0012


def syst_d(x):
    return 0.0009 - x * 0.0009


def _add_e155_points(ds, kin, ww, Qbins, xbins, proc_id, s, thFactor_is_x, syst_fn):
    for i, row in enumerate(kin):
        x_val = row[0]
        point = dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=proc_id,
            s=s, Q_min=Qbins[i][0], Q_max=Qbins[i][1], Q_avg=sqrt(row[1]),
            x_min=xbins[i][0], x_max=xbins[i][1], x_avg=x_val,
            thFactor=x_val if thFactor_is_x else 1.,
            xSec=row[2] - x_val * ww[i][0],
            uncorrErr_0=(row[3] - row[4]) / 2.,
        )
        if syst_fn is not None:
            point["uncorrErr_1"] = syst_fn(x_val)
            point["uncorrErr_2"] = x_val * ww[i][1]
        else:
            point["uncorrErr_1"] = x_val * ww[i][1]
        ds.add_point(point)


#%% -- E155-29.p (2003, 29.1 GeV)
kin1 = _read_csv_rows(path_to_data + "E155/2003/Table1.csv", 29, 8)
kin2 = _read_csv_rows(path_to_data + "E155/2003/Table2.csv", 28, 7)
kin3 = _read_csv_rows(path_to_data + "E155/2003/Table3.csv", 26, 5)
kin = kin1 + kin2 + kin3
ww1 = _read_csv_rows(path_to_ww + "E155_2003_1P(29.1GeV).csv", 0, 8)
ww2 = _read_csv_rows(path_to_ww + "E155_2003_2P(29.1GeV).csv", 0, 7)
ww3 = _read_csv_rows(path_to_ww + "E155_2003_3P(29.1GeV).csv", 0, 5)
ww_29p = ww1 + ww2 + ww3
Qbins_29p = _merge_bins([kin1, kin2, kin3], 1, _guess_bin_sqrt)
xbins_29p = _merge_bins([kin1, kin2, kin3], 0, _guess_bin_linear)
s = 2 * 29.1 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="E155-29.p", comment="Taken from tables 1,2, 3. Bin in Q and x guessed by us.",
                    reference="hep-ex/0204028v1", normErr=[], isNormalized=False)
_add_e155_points(ds, kin, ww_29p, Qbins_29p, xbins_29p, 100, s, thFactor_is_x=True, syst_fn=syst_p)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E155-29.d (2003, 29.1 GeV)
kin1 = _read_csv_rows(path_to_data + "E155/2003/Table1.csv", 56, 8)
kin2 = _read_csv_rows(path_to_data + "E155/2003/Table2.csv", 53, 7)
kin3 = _read_csv_rows(path_to_data + "E155/2003/Table3.csv", 47, 5)
kin = kin1 + kin2 + kin3
ww1 = _read_csv_rows(path_to_ww + "E155_2003_1D(29.1GeV).csv", 0, 8)
ww2 = _read_csv_rows(path_to_ww + "E155_2003_2D(29.1GeV).csv", 0, 7)
ww3 = _read_csv_rows(path_to_ww + "E155_2003_3D(29.1GeV).csv", 0, 5)
ww_29d = ww1 + ww2 + ww3
Qbins_29d = _merge_bins([kin1, kin2, kin3], 1, _guess_bin_sqrt)
xbins_29d = _merge_bins([kin1, kin2, kin3], 0, _guess_bin_linear)
s = 2 * 29.1 * M_deuteron + M_deuteron**2

ds = DataSet.empty("G2", name="E155-29.d", comment="Taken from tables 1,2, 3. Bin in Q and x guessed by us.",
                    reference="hep-ex/0204028v1", normErr=[], isNormalized=False)
_add_e155_points(ds, kin, ww_29d, Qbins_29d, xbins_29d, 102, s, thFactor_is_x=True, syst_fn=syst_d)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E155-32.p (2003, 32.3 GeV; reuses E155-29.p's WW terms, per old parsing's own comment)
kin1 = _read_csv_rows(path_to_data + "E155/2003/Table4.csv", 29, 8)
kin2 = _read_csv_rows(path_to_data + "E155/2003/Table5.csv", 28, 7)
kin3 = _read_csv_rows(path_to_data + "E155/2003/Table6.csv", 26, 5)
kin = kin1 + kin2 + kin3
Qbins_32p = _merge_bins([kin1, kin2, kin3], 1, _guess_bin_sqrt)
xbins_32p = _merge_bins([kin1, kin2, kin3], 0, _guess_bin_linear)
s = 2 * 32.3 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="E155-32.p", comment="Taken from tables 4, 5, 6. Bin in Q and x guessed by us.",
                    reference="hep-ex/0204028v1", normErr=[], isNormalized=False)
_add_e155_points(ds, kin, ww_29p, Qbins_32p, xbins_32p, 100, s, thFactor_is_x=True, syst_fn=syst_p)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E155-32.d (2003, 32.3 GeV; reuses E155-29.d's WW terms)
kin1 = _read_csv_rows(path_to_data + "E155/2003/Table4.csv", 56, 8)
kin2 = _read_csv_rows(path_to_data + "E155/2003/Table5.csv", 53, 7)
kin3 = _read_csv_rows(path_to_data + "E155/2003/Table6.csv", 47, 5)
kin = kin1 + kin2 + kin3
Qbins_32d = _merge_bins([kin1, kin2, kin3], 1, _guess_bin_sqrt)
xbins_32d = _merge_bins([kin1, kin2, kin3], 0, _guess_bin_linear)
s = 2 * 32.3 * M_deuteron + M_deuteron**2

### old comment says "tables 1,2, 3" (copy-paste leftover from the 29.1 GeV deuteron block), ported verbatim
ds = DataSet.empty("G2", name="E155-32.d", comment="Taken from tables 1,2, 3. Bin in Q and x guessed by us.",
                    reference="hep-ex/0204028v1", normErr=[], isNormalized=False)
_add_e155_points(ds, kin, ww_29d, Qbins_32d, xbins_32d, 102, s, thFactor_is_x=False, syst_fn=syst_d)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E155-38.p (1999, 38.8 GeV, no syst. uncertainty)
kin1 = _read_csv_rows(path_to_data + "E155/1999/Table1.csv", 30, 10)
kin2 = _read_csv_rows(path_to_data + "E155/1999/Table2.csv", 27, 7)
kin3 = _read_csv_rows(path_to_data + "E155/1999/Table3.csv", 25, 5)
kin = kin1 + kin2 + kin3
ww1 = _read_csv_rows(path_to_ww + "E155_1999_1P.csv", 0, 10)
ww2 = _read_csv_rows(path_to_ww + "E155_1999_2P.csv", 0, 7)
ww3 = _read_csv_rows(path_to_ww + "E155_1999_3P.csv", 0, 5)
ww_38p = ww1 + ww2 + ww3
Qbins_38p = _merge_bins([kin1, kin2, kin3], 1, _guess_bin_sqrt)
xbins_38p = _merge_bins([kin1, kin2, kin3], 0, _guess_bin_linear)
s = 2 * 38.8 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="E155-38.p", comment="Taken from tables 1, 2, 3. Bin in Q and x guessed by us. No syst. uncert.",
                    reference="hep-ex/9901006v1", normErr=[], isNormalized=False)
_add_e155_points(ds, kin, ww_38p, Qbins_38p, xbins_38p, 100, s, thFactor_is_x=True, syst_fn=None)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E155-38.d (1999, 38.8 GeV, no syst. uncertainty)
kin1 = _read_csv_rows(path_to_data + "E155/1999/Table1.csv", 62, 10)
kin2 = _read_csv_rows(path_to_data + "E155/1999/Table2.csv", 53, 7)
kin3 = _read_csv_rows(path_to_data + "E155/1999/Table3.csv", 47, 5)
kin = kin1 + kin2 + kin3
ww1 = _read_csv_rows(path_to_ww + "E155_1999_1D.csv", 0, 10)
ww2 = _read_csv_rows(path_to_ww + "E155_1999_2D.csv", 0, 7)
ww3 = _read_csv_rows(path_to_ww + "E155_1999_3D.csv", 0, 5)
ww_38d = ww1 + ww2 + ww3
Qbins_38d = _merge_bins([kin1, kin2, kin3], 1, _guess_bin_sqrt)
xbins_38d = _merge_bins([kin1, kin2, kin3], 0, _guess_bin_linear)
s = 2 * 38.8 * M_deuteron + M_deuteron**2

### thFactor=x here too (unlike E155-29.d/E155-32.d, which use thFactor=1) -- ported verbatim, see module NOTE
ds = DataSet.empty("G2", name="E155-38.d", comment="Taken from tables 1, 2, 3. Bin in Q and x guessed by us. No syst. uncert.",
                    reference="hep-ex/9901006v1", normErr=[], isNormalized=False)
_add_e155_points(ds, kin, ww_38d, Qbins_38d, xbins_38d, 102, s, thFactor_is_x=True, syst_fn=None)
ds.save_csv(path_to_save + ds.name + ".csv")
