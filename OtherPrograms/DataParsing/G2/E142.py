"""
Parsing of E142 g2 (transverse spin structure function) data, neutron.

Source: /data/arTeMiDe_Repository/data/g2Tables/E142/Table4.csv
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/E142_N.csv (WW term)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "E142 neutron" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins,
WW-term subtraction, and error treatment kept identical to the old parsing
(the user confirmed this script's physics/process definitions are already
in the modern standard; this is a pure reformatting port, nothing changed).

NOTE: process code. Old parsing used a single `process` integer (100=p,
101=n, 102=d). Per user instruction, this maps directly onto Cerynia's G2
schema as proc_id (unchanged value), with ps_def=1, h_1=1 fixed (G2 has no
h_2 -- inclusive DIS, no produced hadron).

NOTE: xSec = raw_g2 - WW_value (Wandzura-Wilczek/leading-twist term,
precomputed offline and stored in WW_MAPPDFpol/E142_N.csv, subtracted
directly, not x-scaled -- unlike some other G2 sub-cases in this category).
uncorrErr = [stat(symmetrized), syst(symmetrized), WW_stat]. thFactor=1.

NOTE: Q bin is "guessed" (not given directly in the raw table) via the old
script's GuessBin_Bottom/Middle/Top helpers, ported verbatim below: first
row uses Bottom, last row uses Top, interior rows use Middle, operating on
sqrt(Q^2) using the neighboring rows' Q^2 values.

NOTE: normErr=[] -- the old script computes a `lumUncertainty` value for
every G2 dataset but never actually appends it to any DataSet's normErr
(confirmed: every old G2 CSV has "Number of norm.errors,0"). Ported as
found, not a per-dataset choice.

Dataset (arXiv:hep-ex/9610007), 8 points:
    E142.n
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

ps_def, h_1 = 1, 1


def _guess_bin_sqrt(data, i, col):
    """
    Ported verbatim from old parsing's GuessBin_Bottom/Middle/Top("Q",...):
    guesses [lo,hi] for sqrt(data[i][col]) from neighboring rows' values.
    First row uses "bottom" (no left neighbor), last row uses "top" (no
    right neighbor), interior rows use "middle" (both neighbors).
    """
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
### columns: [0]=x_avg, [1]=x_min, [2]=x_max, [3]=Q2_avg, [4]=g2(raw), [5]=stat+, [6]=stat-, [7]=syst+, [8]=syst-
data = _read_csv_rows(path_to_data + "E142/Table4.csv", 25, 8)
ww = _read_csv_rows(path_to_ww + "E142_N.csv", 0, 8)

s = 2 * 22.367 * M_neutron + M_neutron**2

ds = DataSet.empty("G2", name="E142.n",
                    comment="Taken from table 4. Bin in Q guessed by us. We take the three incident energies presented for the data and assume the values computed are for the average energy from those 3",
                    reference="hep-ex/9610007", normErr=[], isNormalized=False)

for i, row in enumerate(data):
    Q_min, Q_max = _guess_bin_sqrt(data, i, 3)
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=101,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[3]),
        x_min=row[1], x_max=row[2], x_avg=row[0],
        thFactor=1.,
        xSec=row[4] - ww[i][0],
        uncorrErr_0=(row[5] - row[6]) / 2.,
        uncorrErr_1=(row[7] - row[8]) / 2.,
        uncorrErr_2=ww[i][1],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
