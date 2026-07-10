"""
Parsing of HERMES (2020 publication) double-spin asymmetry
A_LT^cos(phi-phis)/sqrt(1-e^2) for pi+/pi-/k+/k-, proton target, 3D
(x,z,pT)-binned.

Source: /data/arTeMiDe_Repository/data/HERMES-SSA/SFA/hermesTMDs_data__{pip,pim,kp,km}_SFA_3D.txt
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/wgt/*.csv

Ported from DataProcessor/OtherPrograms/DataParsing/ReadALT_HERMES.py -- s,
x/z/pT/Q bins, cuts, and error treatment kept identical to the old parsing.
Reads from the SAME raw files as Sivers/hermes3D.py, just a different
asymmetry block within each file (each file holds 8 different transverse/
double-spin asymmetries; this one picks out the "A_LT^cos(phi-phis)/
sqrt(1-e^2)" block by its header text, same 64-point (x,z,pT) 3D binning
as the Sivers block, same Qbounds() kinematic-Q-window helper).

NOTE: process code, already correct in old parsing for pi+/pi-/k+ (matches
the now-familiar [1,1,h_2,13001] primary / [1,1,h_2,2001] weight scheme,
h_1=1 proton -- genuinely correct here, no target-species bug like the
COMPASS cases). h_2: pi+=1, pi-=-1, k+=2.

NOTE: k- bugfix. Old parsing's 'km' branch used h_2=-1 (the SAME code as
pi-), instead of -2 (which is what k+ uses +2 for, and what every other
pi/k hadron-code pairing in this migration uses: pi+-=+-1, k+-=+-2) --
almost certainly a copy-paste error from the 'pim' branch just above it in
the old script. Confirmed and fixed per user instruction: h_2=-2 for k-.

NOTE: naming. Old "hermes3D.ALT.<hadron>" -> "hermes3D.wgt.<hadron>"
(".ALT." -> ".wgt.", matching the category-folder rename already applied
to compass16.wgt.*), not explicitly requested for this sub-case but
applied for naming consistency across the wgt category -- flag if a
different convention is wanted.

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention). Real
per-point x/z/pT bin edges (unlike hermes09-style fixed shared bins),
same as Sivers/hermes3D.py.

Datasets, 4 total (reference 2007.07755, 64 points each):
    hermes3D.wgt.pi+, hermes3D.wgt.pi-, hermes3D.wgt.k+, hermes3D.wgt.k-
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/HERMES-SSA/SFA/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/wgt/"

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

ps_def, h_1 = 1, 1
proc_id, proc_id_weight = 13001, 2001

START_MARKER = "A_LT^cos(phi-phis)/sqrt(1-e^2) SFA DSA"
END_MARKER = "in 3D bins of x, z, Pt"


def Qbounds(xMin, xMax):
    """
    Kinematic Q bounds given an x window (ported verbatim from old parsing,
    same as Sivers/hermes3D.py):
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


def _read_hermes3d_alt_table(path):
    """
    Find the "A_LT^cos(phi-phis)/sqrt(1-e^2)" block (by its unique header
    text) and read its 64 (x,z,pT)-binned rows, stopping at the next
    block's header (also containing "in 3D bins of x, z, Pt").
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    reading = False
    for line in lines:
        if START_MARKER in line:
            reading = True
            continue
        if reading:
            if END_MARKER in line:
                break
            if "#" in line or line == "":
                continue
            parts = line.split()
            x_bin = [float(v) for v in parts[0].split("<x<")]
            z_bin = [float(v) for v in parts[1].split("<z<")]
            pT_bin = [float(v) for v in parts[2].split("<Pt<")]
            rest = [float(v) for v in parts[3:12]]
            rows.append((x_bin, z_bin, pT_bin, rest))
    return rows


def _add_hermes3d_alt_points(ds, rows, h_2, m_product):
    for i, (x_bin, z_bin, pT_bin, rest) in enumerate(rows):
        ### rest = [Q2_avg, x_avg, y_avg(unused), z_avg, pT_avg, e_avg(unused), xSec, stat, syst]
        Q_min, Q_max = Qbounds(x_bin[0], x_bin[1])
        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=h_2, proc_id_weight=proc_id_weight,
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
rows = _read_hermes3d_alt_table(path_to_data + "hermesTMDs_data__pip_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.wgt.pi+",
                    comment="HERMES 3D DSA-A_LT pi+ from p. The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_alt_points(ds, rows, 1, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-
rows = _read_hermes3d_alt_table(path_to_data + "hermesTMDs_data__pim_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.wgt.pi-",
                    comment="HERMES 3D DSA-A_LT pi- from p. The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_alt_points(ds, rows, -1, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+
rows = _read_hermes3d_alt_table(path_to_data + "hermesTMDs_data__kp_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.wgt.k+",
                    comment="HERMES 3D DSA-A_LT k+ from p. The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_alt_points(ds, rows, 2, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-
### fixed vs old parsing: h_2=-2 (old code had -1, a copy-paste bug from the pi- branch), see module NOTE
rows = _read_hermes3d_alt_table(path_to_data + "hermesTMDs_data__km_SFA_3D.txt")
ds = DataSet.empty("SIDIS", name="hermes3D.wgt.k-",
                    comment="HERMES 3D DSA-A_LT k- from p. The data MUST be evaluated at a point",
                    reference="2007.07755", normErr=normErr, isNormalized=False)
_add_hermes3d_alt_points(ds, rows, -2, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")
