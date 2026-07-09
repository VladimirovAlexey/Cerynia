"""
Parsing of E772 (fixed-target, mu+mu- Drell-Yan) invariant cross-section data.

Source: /data/arTeMiDe_Repository/data/E772/E772_4to9.dat, E772_11to15.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/E772.csv

Ported from the "E772" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q-range, y-range, thFactor, error split) are kept identical
to the old parsing. Both raw files feed into the single "E772" dataset, as in
old parsing.

NOTE: thFactor is the fixed-target "invariant cross section" Jacobian,
1/(qT_max**2-qT_min**2)/(y_max-y_min)*ffactor -- structurally different from
the collider bin-integrated k/(qT_max-qT_min) form. Per user confirmation,
atmdeFactor still applies here: the old expression is multiplied by
atmdeFactor=(Q_max-Q_min)*(y_max-y_min)*(qT_max-qT_min), same as every other
DY case.

CORRECTED (per audit of DataProcessor.git history): the originally-ported
ReadDYdataFiles.py block used a plain constant 0.3183098861838 (=1/pi) in
place of ffactor -- missing a Feynman-x Jacobian. Commit "Corrected thFactor
in E605 and E772" (2026-01-05, DataProcessor's own repo) fixed this in a
parallel script (ReadDYdataFiles(low-energy).py) that never got merged back
into ReadDYdataFiles.py. ffactor = sqrt(4*(Q_avg**2+qT_avg**2)/s + xF**2) *
0.3183098861838, where xF is this dataset's Feynman-x (here the midpoint of
the y-window, (y_min+y_max)/2=0.2 -- the fix's own variable was misleadingly
named "<y>", but the accompanying comment makes clear it means xF, per
xF = x1+x2 = sqrt(4*(Q^2+qT^2)/s + xF^2) at fixed-target kinematics).

Also: this raw table's x/xlow/xhigh columns are degenerate (all equal to the
bin center -- no real bin-width info), so qT bounds are manually injected as
center +- 0.125 GeV, exactly as old parsing did.

Dataset:
    E772 -- Phys.Rev.D 50 (1994) 3-38 + Erratum D60 (1999) 119903
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"


def _read_hepdata_table(path, n_header):
    """
    Read a tab-separated HEPData .dat table.
    Each data row is "\\t x xlow xhigh y dy+ dy- [dy+ dy- ...]" -- the leading
    tab produces an empty first field, which is dropped. Returns a list of
    float rows [x, xlow, xhigh, y, dy+, dy-, ...].
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    lines = lines[n_header:-1]  # drop header block and trailing blank line
    return [[float(v) for v in line.split("\t")[1:]] for line in lines]


#%%
# ============================================================================
# E772 -- Phys.Rev.D 50 (1994) 3-38 + Erratum D60 (1999) 119903
# ============================================================================
### 800 GeV fixed-target beam momentum -> sqrt(s)
s = 38.76**2
### Process is virtual-photon DY for p-deuteron (target treated as isoscalar proton)
## ps_def, h_1, h_2, proc_id = 2, 1, 1, 102 ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 2, 1, 12, 1 ##(>6.04 definition)
### XL : 0.1 TO 0.3, stated in both file headers (used as a y-window in old parsing)
y_min, y_max = 0.1, 0.3
### 10% normalization uncertainty (as in old parsing)
normErr = [0.10]
### Feynman-x for this dataset = midpoint of the y-window (as in old parsing's fix)
xF = (y_min + y_max) * 0.5

ds = DataSet.empty("DY", name="E772", comment="E772 data",
                    reference="Phys.Rev.D 50 (1994) 3-38 + Erratum D60 (1999) 119903",
                    normErr=normErr, isNormalized=False)

### File 1: mass sub-bins M=5.5,6.5,7.5,8.5 GeV -> Q-windows [5,6],[6,7],[7,8],[8,9]
rows = _read_hepdata_table(path_to_data + "E772/E772_4to9.dat", n_header=8)

for j in range(4):
    Q_min, Q_max = float(5 + j), float(6 + j)
    for i, row in enumerate(rows):
        x = row[0]
        qT_min, qT_max = x - 0.125, x + 0.125  ### x/xlow/xhigh are degenerate; bin width injected manually
        xSec = row[3 + 3 * j]
        dyp, dym = row[4 + 3 * j], row[5 + 3 * j]
        if xSec == -50:  # missing-data sentinel (as in old parsing)
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ### invariant cross-section Jacobian, corrected for Feynman-x (as in the
        ### 2026-01-05 fix); 0.3183098861838 = 1/pi
        Q_avg = (Q_min + Q_max) * 0.5
        qT_avg = (qT_min + qT_max) * 0.5
        ffactor = sqrt(4. * (Q_avg**2 + qT_avg**2) / s + xF**2) * 0.3183098861838

        ds.add_point(dict(
            id=f"E772.{int(Q_min)}Q{int(Q_max)}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * ffactor,
            includeCuts=False,
            xSec=xSec,
            uncorrErr_0=(dyp - dym) / 2.,
        ))

### File 2: mass sub-bins M=11.5,12.5,13.5,14.5 GeV -> Q-windows [11,12],[12,13],[13,14],[14,15]
rows = _read_hepdata_table(path_to_data + "E772/E772_11to15.dat", n_header=8)

for j in range(4):
    Q_min, Q_max = float(11 + j), float(12 + j)
    for i, row in enumerate(rows):
        x = row[0]
        qT_min, qT_max = x - 0.125, x + 0.125
        xSec = row[3 + 3 * j]
        dyp, dym = row[4 + 3 * j], row[5 + 3 * j]
        if xSec == -50:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ### invariant cross-section Jacobian, corrected for Feynman-x (as in the
        ### 2026-01-05 fix); 0.3183098861838 = 1/pi
        Q_avg = (Q_min + Q_max) * 0.5
        qT_avg = (qT_min + qT_max) * 0.5
        ffactor = sqrt(4. * (Q_avg**2 + qT_avg**2) / s + xF**2) * 0.3183098861838

        ds.add_point(dict(
            id=f"E772.{int(Q_min)}Q{int(Q_max)}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * ffactor,
            includeCuts=False,
            xSec=xSec,
            uncorrErr_0=(dyp - dym) / 2.,
        ))

ds.save_csv(path_to_save + ds.name + ".csv")
