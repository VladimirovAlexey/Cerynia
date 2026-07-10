"""
Parsing of HERMES (2020 publication) transverse single-spin asymmetry
A_UT^Sivers for pi+/pi-/k+/k-, proton target, 3D (x,z,pT)-binned.

Source: /data/arTeMiDe_Repository/data/HERMES-SSA/SFA/hermesTMDs_data__{pip,pim,kp,km}_SFA_3D.txt
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from the "hermes.sivers.<hadron>.3d" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSSA_HERMES.py -- s, x/z/pT/Q
bins, cuts, and error treatment kept identical to the old parsing. This is
a separate, newer measurement from the 2009-publication Q-integrated/
Q<2/Q>2-sliced data -- see hermes09.py for those (and for the Qbounds()
kinematic-Q-window helper, reused here unchanged).

NOTE: process code. Confirmed correct as-is by the user -- primary process
[1,1,h_2,12001] and weight process [1,1,h_2,2001] (h_2: pi+=1, pi-=-1,
k+=2, k-=-2), matching the old script's "process"/"weightProcess" fields
exactly. h_1=1 (proton) is correct (HERMES uses a proton target).

NOTE: naming convention, given by the user: "hermes.sivers.<hadron>.3d" ->
"hermes3D.sivers.<hadron>" (dropping the ".3d" suffix, moving the "3D"
marker into the name prefix), matching the naming style already used for
hermes3D.* in the SIDIS category (see feedback_sidis_conventions).

NOTE: unlike hermes09.py, x/z/pT bins here are REAL per-point bin edges
(read directly from the raw table's "lo<x<hi" fields), not the fixed
shared bin-edge lists used for the 2009 data -- still evaluated at a point
per the old dataset comment, but the bin itself is meaningful (used for
the Qbounds() Q-window calculation).

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention).

Datasets, 4 total:
    hermes3D.sivers.pi+, hermes3D.sivers.pi-, hermes3D.sivers.k+, hermes3D.sivers.k-
    64 points each, reference 2007.07755
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/HERMES-SSA/SFA/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_proton = 0.938
m_pion = 0.139
m_kaon = 0.494

### 27.6 GeV beam on a fixed proton target (as in old parsing)
s = 2 * 27.6 * 0.938 + 0.938**2
### 7.3% overall scale uncertainty (as in old parsing)
normErr = [0.073]
includeCuts = True
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.95, 10., 10000.]


def Qbounds(xMin, xMax):
    """
    Kinematic Q bounds given an x window (ported verbatim from old parsing,
    same as hermes09.py):
    Qmin^2 = MAX(Q2min, xMin*yMin*(s-M^2), xMin/(1-xMin)*(W2min-M^2))
    Qmax^2 = MIN(Q2max, xMax*yMax*(s-M^2), xMax/(1-xMax)*(W2max-M^2))
    """
    Q2min, Q2max = 1., 10000.
    WM2min, WM2max = 10. - 0.938**2, 10000. - 0.938**2
    yMin, yMax = 0.1, 0.95
    sM2 = 2 * 27.6 * 0.938
    if xMax < 1:
        return [sqrt(max(Q2min, xMin * yMin * sM2, xMin / (1 - xMin) * WM2min)),
                sqrt(min(Q2max, xMax * yMax * sM2, xMax / (1 - xMax) * WM2max))]
    else:
        return [sqrt(max(Q2min, xMin * yMin * sM2, xMin / (1 - xMin) * WM2min)),
                sqrt(min(Q2max, xMax * yMax * sM2, 1000 * WM2max))]


def _read_hermes3d_table(path):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    for line in lines[86:150]:
        parts = line.split()
        x_bin = [float(v) for v in parts[0].split("<x<")]
        z_bin = [float(v) for v in parts[1].split("<z<")]
        pT_bin = [float(v) for v in parts[2].split("<Pt<")]
        ### rest = [Q2_avg, x_avg, y_avg(unused), z_avg, pT_avg, e_avg(unused), xSec, stat, syst]
        rest = [float(v) for v in parts[3:12]]
        rows.append((x_bin, z_bin, pT_bin, rest))
    return rows


def _add_hermes3d_points(ds, rows, h_2, m_product):
    for i, (x_bin, z_bin, pT_bin, rest) in enumerate(rows):
        Q_min, Q_max = Qbounds(x_bin[0], x_bin[1])
        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=1, h_1=1, h_2=h_2, proc_id=12001,
            ps_def_weight=1, h_1_weight=1, h_2_weight=h_2, proc_id_weight=2001,
            s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(rest[0]),
            x_min=x_bin[0], x_max=x_bin[1], x_avg=rest[1],
            z_min=z_bin[0], z_max=z_bin[1], z_avg=rest[3],
            pT_min=pT_bin[0], pT_max=pT_bin[1], pT_avg=rest[4],
            M_target=M_proton, M_product=m_product,
            includeCuts=includeCuts,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            thFactor=1.,
            xSec=rest[6],
            uncorrErr_0=rest[7],
            uncorrErr_1=rest[8],
        ))


#%% -- pi+
rows = _read_hermes3d_table(path_to_data + "hermesTMDs_data__pip_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.sivers.pi+",
                    comment="HERMES SSA-Sivers pi+ (3d-data). The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_points(ds, rows, 1, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-
rows = _read_hermes3d_table(path_to_data + "hermesTMDs_data__pim_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.sivers.pi-",
                    comment="HERMES SSA-Sivers pi- (3d-data). The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_points(ds, rows, -1, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+
rows = _read_hermes3d_table(path_to_data + "hermesTMDs_data__kp_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.sivers.k+",
                    comment="HERMES SSA-Sivers k+ (3d-data). The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_points(ds, rows, 2, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-
rows = _read_hermes3d_table(path_to_data + "hermesTMDs_data__km_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.sivers.k-",
                    comment="HERMES SSA-Sivers k- (3d-data). The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_points(ds, rows, -2, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")
