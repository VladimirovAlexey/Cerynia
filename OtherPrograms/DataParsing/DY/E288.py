"""
Parsing of E288 (fixed-target, mu+mu- Drell-Yan) invariant cross-section data.

Source: /data/arTeMiDe_Repository/data/E288/E288_200.dat, E288_300.dat, E288_400.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "E288" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q-range, y-range, thFactor, error split) are kept identical
to the old parsing.

NOTE: this dataset's old thFactor is the fixed-target "invariant cross
section" Jacobian, 1/(qT_max**2-qT_min**2)/(y_max-y_min)*0.3183098861838*1000,
structurally different from the collider bin-integrated k/(qT_max-qT_min)
form seen in Tevatron/ATLAS/CMS/LHCb. Per user confirmation, atmdeFactor
still applies here -- the old expression (whatever its shape) is multiplied
by atmdeFactor=(Q_max-Q_min)*(y_max-y_min)*(qT_max-qT_min), same as every
other DY case. Also note the old code's inline comment said "0.001 for nb"
but the actual factor applied was *1000 (not *0.001); ported verbatim as
coded (now folded into the atmdeFactor-multiplied form below) -- flagged for
the user to review, not corrected.

CORRECTED (per audit of DataProcessor.git history): ps_def was 2 in the
originally-ported ReadDYdataFiles.py block, but a later, undated-in-that-file
fix (commit "Corrections in E228", 2026-01-02, in DataProcessor's own repo)
changed it to 1 for all three sub-datasets. Applied here to both the ">6.04"
(h_2=63029) and the historical "<6.04" comment line, for consistency.

Datasets (Phys.Rev.D 23 (1981) 604). Old parsing named these "E228-*" --
a typo for "E288-*", fixed here per user instruction:
    E288-200 -- 200 GeV incident beam
    E288-300 -- 300 GeV incident beam
    E288-400 -- 400 GeV incident beam
"""

import sys
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
# E288, 200 GeV beam -- Phys.Rev.D 23 (1981) 604
# ============================================================================
### "(AV empty data points filled by -50)", stated in the file header
rows = _read_hepdata_table(path_to_data + "E288/E288_200.dat", n_header=9)

### 200 GeV fixed-target beam momentum -> sqrt(s)
s = 19.42**2
### Process is virtual-photon DY for p-copper (target treated as isoscalar proton)
##ps_def, h_1, h_2, proc_id = 1, 1, 1, 101  ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 1, 1, 63029, 1  ##(>6.04 definition)

### YRAP(P=3 4,RF=CM) : 0.4, stated in the file header (fixed rapidity bin, +-0.3 half-width as in old parsing)
y_min, y_max = 0.1, 0.7
### 25% normalization uncertainty (as in old parsing)
normErr = [0.25]

ds = DataSet.empty("DY", name="E288-200", comment="E288 (200) data", reference="Phys.Rev.D 23 (1981) 604",
                    normErr=normErr, isNormalized=False)

### 7 mass sub-bins, W=4-5,5-6,...,10-11 GeV, one xSec+err column-block per bin
### (reduced offsets 3+3*j, 4+3*j, 5+3*j)
for j in range(7):
    Q_min, Q_max = float(4 + j), float(5 + j)
    for i, row in enumerate(rows):
        qT_min, qT_max = row[1], row[2]
        xSec = row[3 + 3 * j]
        dyp, dym = row[4 + 3 * j], row[5 + 3 * j]
        if xSec == -50:  # missing-data sentinel (as in old parsing)
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"E288-200.{int(Q_min)}Q{int(Q_max)}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            ### invariant cross-section Jacobian (as in old parsing); 0.3183098861838 = 1/pi
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * 0.3183098861838 * 1000,
            includeCuts=False,
            xSec=xSec,
            uncorrErr_0=(dyp - dym) / 2.,
        ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# E288, 300 GeV beam -- Phys.Rev.D 23 (1981) 604
# ============================================================================
rows = _read_hepdata_table(path_to_data + "E288/E288_300.dat", n_header=9)

### 300 GeV fixed-target beam momentum -> sqrt(s)
s = 23.73**2
## ps_def, h_1, h_2, proc_id = 1, 1, 1, 101 ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 1, 1, 63029, 1  ##(>6.04 definition)

### YRAP(P=3 4,RF=CM) : 0.21, stated in the file header, +-0.3 half-width as in old parsing
y_min, y_max = 0.21 - 0.3, 0.21 + 0.3
normErr = [0.25]

ds = DataSet.empty("DY", name="E288-300", comment="E288 (300) data", reference="Phys.Rev.D 23 (1981) 604",
                    normErr=normErr, isNormalized=False)

### 8 mass sub-bins, W=4-5,5-6,...,11-12 GeV
for j in range(8):
    Q_min, Q_max = float(4 + j), float(5 + j)
    for i, row in enumerate(rows):
        qT_min, qT_max = row[1], row[2]
        xSec = row[3 + 3 * j]
        dyp, dym = row[4 + 3 * j], row[5 + 3 * j]
        if xSec == -50:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"E288-300.{int(Q_min)}Q{int(Q_max)}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * 0.3183098861838 * 1000,
            includeCuts=False,
            xSec=xSec,
            uncorrErr_0=(dyp - dym) / 2.,
        ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# E288, 400 GeV beam -- Phys.Rev.D 23 (1981) 604
# ============================================================================
rows = _read_hepdata_table(path_to_data + "E288/E288_400.dat", n_header=9)

### 400 GeV fixed-target beam momentum -> sqrt(s)
s = 27.43**2
## ps_def, h_1, h_2, proc_id = 1, 1, 1, 101 ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 1, 1, 63029, 1  ##(>6.04 definition)

### YRAP(P=3 4,RF=CM) : 0.03, stated in the file header, +-0.3 half-width as in old parsing
y_min, y_max = 0.03 - 0.3, 0.03 + 0.3
normErr = [0.25]

ds = DataSet.empty("DY", name="E288-400", comment="E288 (400) data", reference="Phys.Rev.D 23 (1981) 604",
                    normErr=normErr, isNormalized=False)

### 9 mass sub-bins, W=5-6,6-7,...,13-14 GeV
for j in range(9):
    Q_min, Q_max = float(5 + j), float(6 + j)
    for i, row in enumerate(rows):
        qT_min, qT_max = row[1], row[2]
        xSec = row[3 + 3 * j]
        dyp, dym = row[4 + 3 * j], row[5 + 3 * j]
        if xSec == -50:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"E288-400.{int(Q_min)}Q{int(Q_max)}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * 0.3183098861838 * 1000,
            includeCuts=False,
            xSec=xSec,
            uncorrErr_0=(dyp - dym) / 2.,
        ))

ds.save_csv(path_to_save + ds.name + ".csv")
