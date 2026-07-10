"""
Parsing of JLab (6 GeV, Hall A, polarized 3He / effective-neutron target)
double-spin asymmetry A_LT for pi+/pi-.

Source: hardcoded literal values in the old parsing script (see NOTE below) --
        no raw data file is read; values are taken from a table printed at
        the end of arXiv:1108.0489.
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/wgt/*.csv

Ported from DataProcessor/OtherPrograms/DataParsing/ReadALT_JLab.py -- s,
x/z/pT/Q bins, cuts, and error treatment kept identical to the old parsing.

NOTE: target is neutron, confirmed by the user (old code's own comment says
"off neutron" for both pi+/pi-). Per user instruction, h_1=11 (neutron),
matching the fix already applied to Sivers/JLab.py -- same experiment
family (JLab Hall A, polarized 3He), same target-code correction.

NOTE: process code. Old parsing used [1,1,h_2,13003] (primary) /
[1,1,h_2,2003] (weight) -- non-uniform proc_id, same pattern as the old
JLab Sivers case (which used 12003/2003 before being corrected to the
standard 12001/2001). Per user instruction, updated to primary
[1,11,h_2,13001] / weight [1,11,h_2,2001] (h_2: 1=pi+, -1=pi-), matching
the uniform wgt proc_id already established by compass16.wgt.*/
hermes3D.wgt.* (13001/2001).

NOTE: naming. Old "JLab6.ALT.<hadron>" -> "JLab6.wgt.<hadron>"
(".ALT."->".wgt."), matching the category-wide rename applied throughout
this migration.

NOTE: all bins except x are FIXED literals (Q=[sqrt(1.4),sqrt(2.7)],
z=[0.5,0.6], pT=[0.24,0.44]), taken directly from the paper text, not
per-point measured values -- ported verbatim. x bins are also fixed,
hand-computed as [0.15+0.05*(n-1), 0.15+0.05*n] for n=1..4 (i.e.
[0.15,0.20], [0.20,0.25], [0.25,0.30], [0.30,0.35]) -- same "data must be
evaluated at a point" situation as Sivers/JLab.py (rough bins, not real
per-point edges).

NOTE: M_target stays the proton mass (0.938) even for this neutron target,
same convention as Sivers/JLab.py (proton/neutron mass differ by ~0.06%,
not corrected here, ported verbatim). normErr=[0.028] ("delution factor
according to polarization uncertainty" per old code's comment) -- ported
verbatim.

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention).

Datasets, 2 total (reference 1108.0489, 4 points each):
    JLab6.wgt.pi+, JLab6.wgt.pi-
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/wgt/"

M_proton = 0.938
m_pion = 0.139

### 5.9 GeV electron beam on a fixed nucleon target (as in old parsing, same as Sivers/JLab.py)
s = 2 * 5.9 * 0.938 + 0.938**2
### 2.8% dilution/polarization-uncertainty factor (as in old parsing)
normErr = [0.028]
includeCuts = False
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.9, 10., 10000.]

### neutron target (confirmed by the user, matching Sivers/JLab.py)
ps_def, h_1 = 1, 11
proc_id, proc_id_weight = 13001, 2001

### fixed bins from the paper text (as in old parsing)
Q_min, Q_max = sqrt(1.4), sqrt(2.7)
z_min, z_max = 0.5, 0.6
pT_min, pT_max = 0.24, 0.44

### table from arXiv:1108.0489, [x, Q^2, z, pT, A_pi+, un1, un2, A_pi-, un1, un2]
TABLE = [
    [0.156, 1.38, 0.50, 0.43,  0.02, 0.11, 0.03, 0.10, 0.07, 0.03],
    [0.206, 1.76, 0.52, 0.38,  0.04, 0.13, 0.06, 0.18, 0.11, 0.06],
    [0.265, 2.16, 0.54, 0.32, -0.13, 0.11, 0.05, 0.10, 0.07, 0.02],
    [0.349, 2.68, 0.58, 0.24, -0.27, 0.18, 0.13, 0.18, 0.10, 0.05],
]


def _add_jlab6_points(ds, h_2, xSec_col, err_cols):
    for i, row in enumerate(TABLE):
        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=h_2, proc_id_weight=proc_id_weight,
            s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[1]),
            x_min=0.15 + 0.05 * i, x_max=0.15 + 0.05 * (i + 1), x_avg=row[0],
            z_min=z_min, z_max=z_max, z_avg=row[2],
            pT_min=pT_min, pT_max=pT_max, pT_avg=row[3],
            M_target=M_proton, M_product=m_pion,
            includeCuts=includeCuts,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            thFactor=1.,
            xSec=row[xSec_col],
            uncorrErr_0=row[err_cols[0]],
            uncorrErr_1=row[err_cols[1]],
        ))


#%% -- pi+
ds = DataSet.empty("SIDIS", name="JLab6.wgt.pi+", comment="JLab at 6 GeV data for A_LT pi+ off neutron",
                    reference="1108.0489", normErr=normErr, isNormalized=False)
_add_jlab6_points(ds, h_2=1, xSec_col=4, err_cols=(5, 6))
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-
ds = DataSet.empty("SIDIS", name="JLab6.wgt.pi-", comment="JLab at 6 GeV data for A_LT pi- off neutron",
                    reference="1108.0489", normErr=normErr, isNormalized=False)
_add_jlab6_points(ds, h_2=-1, xSec_col=7, err_cols=(8, 9))
ds.save_csv(path_to_save + ds.name + ".csv")
