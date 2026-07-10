"""
Parsing of RSS g2 (transverse spin structure function) data, proton.

Source: /data/arTeMiDe_Repository/data/g2Tables/RSSC/2007Exp_P.csv
        /data/arTeMiDe_Repository/data/g2Tables/WW_MAPPDFpol/RSSC.csv (WW term)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/G2/*.csv

Ported from the "RSS" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingG2.py -- s, Q/x bins,
WW-term subtraction, and error treatment kept identical to the old parsing
(the user confirmed this script's physics/process definitions are already
in the modern standard; this is a pure reformatting port, nothing changed).

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction); proc_id
carried over unchanged from the old single `process` integer (100=proton).

NOTE: raw data was hand-digitized from a figure ("Plot Digitalizer", per
old comment, no HEPData table exists) -- 26 rows, no header, columns
[x_avg, g2(raw), error(single, not a +/- pair -- used directly, not
symmetrized)]. x bin is guessed (via _guess_bin_linear, same Bottom/
Middle/Top logic as G2/E142.py/E155.py); Q is a FIXED literal window for
every point (Q_avg=sqrt(1.3), Q=[sqrt(0.8),sqrt(1.4)]) -- old comment says
"Q^2 1.36 GeV^2" but the code uses 1.3, a verbatim discrepancy, not
corrected.

NOTE: xSec = raw_g2 - WW_value (unscaled, same pattern as E142/E154).
uncorrErr = [raw_error(unsymmetrized), WW_stat]. thFactor=1.

NOTE: normErr=[] (matches every other G2 dataset).

Dataset (arXiv:nucl-ex/0608003v3), 26 points:
    RSS.p
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


#%% -- proton
### columns: [0]=x_avg, [1]=g2(raw), [2]=error (single, unsymmetrized)
kin = _read_csv_rows(path_to_data + "RSSC/2007Exp_P.csv", 0, 26)
ww = _read_csv_rows(path_to_ww + "RSSC.csv", 0, 26)
xbins = [_guess_bin_linear(kin, i, 0) for i in range(len(kin))]

### fixed Q window (as in old parsing; see module NOTE on the Q^2=1.3 vs 1.36 discrepancy)
Q_min, Q_max, Q_avg = sqrt(0.8), sqrt(1.4), sqrt(1.3)
s = 2 * 5.755 * M_proton + M_proton**2

ds = DataSet.empty("G2", name="RSS.p",
                    comment="Extracted from FIG 4 (in the paper) with the aid of ''Plot Digitalizer'' because there are no Datasets. Bin in x guessed by us. All measurements at the same energy of Q^2 1.36 GeV^2, we assume there's no bin on Q then.",
                    reference="nucl-ex/0608003v3", normErr=[], isNormalized=False)

for i, row in enumerate(kin):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, proc_id=100,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg,
        x_min=xbins[i][0], x_max=xbins[i][1], x_avg=row[0],
        thFactor=1.,
        xSec=row[1] - ww[i][0],
        uncorrErr_0=row[2],
        uncorrErr_1=ww[i][1],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
