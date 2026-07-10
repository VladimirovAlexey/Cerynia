"""
Parsing of COMPASS (2022 data, ~70%) transverse single-spin asymmetry
A_UT^Sivers for isoscalar h+/h-, deuteron target.

Source: /data/arTeMiDe_Repository/data/COMPASS/2401.00309/Sivers.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from the "compass23.sivers.*" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSSA_Compass23.py -- s, x/z/pT/Q
bins, cuts, and error treatment kept identical to the old parsing.

NOTE: process code. Old parsing used [1,1,h_2,12103] (primary) /
[1,1,h_2,2103] (weight), with h_1=1 (proton) despite the raw file's own
title "Sivers asymmetry (COMPASS deuteron 2022 data (~70%))". Per user
confirmation, this is deuteron data: updated to primary [1,12,h_2,12001] /
weight [1,12,h_2,2001] (h_1=12 deuteron; h_2=12 for h+, -12 for h-, same
isoscalar-hadron code as compass.d.h+/h- in the unpolarized-SIDIS category,
see feedback_sidis_conventions). proc_id=12001/weight=2001 also matches the
now-uniform SIDIS-Sivers process code seen in JLab.py and hermes09.py/
hermes3D.py -- proc_id=12001 is the same across every SIDIS Sivers
sub-case regardless of target/experiment; only h_1/h_2 vary.

NOTE: M_target stays the proton mass (0.938) even for this deuteron target
-- ported verbatim, matching the established per-nucleon SIDIS convention
(see feedback_sidis_conventions: deuteron-target SIDIS structure functions
use the single-nucleon mass). M_product uses the pion mass as a stand-in
for the unidentified isoscalar hadron (same convention as compass.d.h+/h-).

NOTE: only the TOTAL error column (er_tot) is used as uncorrErr_0; the
separate stat error column (er_stat) is read but not stored, matching old
parsing exactly (only one uncorrErr appended).

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention).

Datasets, 6 total (arXiv:2401.00309):
    compass23.sivers.h+.dx, .dz, .dpt
    compass23.sivers.h-.dx, .dz, .dpt
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/COMPASS/2401.00309/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_proton = 0.938
m_pion = 0.139

### 160 GeV beam on a fixed deuteron target (as in old parsing)
s = 2. * 160. * 0.938 + 0.938**2
includeCuts = True
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.9, 25., 10000.]

### deuteron target (per user instruction); h_2: 12=h+, -12=h- (isoscalar hadron)
ps_def, h_1 = 1, 12
proc_id_weight = 2001


def Qbounds(xMin, xMax):
    """
    Kinematic Q bounds given an x window (ported verbatim from old parsing):
    Qmin^2 = MAX(Q2min, xMin*yMin*(s-M^2), xMin/(1-xMin)*(W2min-M^2))
    Qmax^2 = MIN(Q2max, xMax*yMax*(s-M^2), xMax/(1-xMax)*(W2max-M^2))
    """
    Q2min, Q2max = 1., 10000.
    WM2min, WM2max = 25. - 0.938**2, 10000. - 0.938**2
    yMin, yMax = 0.1, 0.9
    sM2 = 2. * 160. * 0.938
    if xMax < 1:
        return [sqrt(max(Q2min, xMin * yMin * sM2, xMin / (1 - xMin) * WM2min)),
                sqrt(min(Q2max, xMax * yMax * sM2, xMax / (1 - xMax) * WM2max))]
    else:
        return [sqrt(max(Q2min, xMin * yMin * sM2, xMin / (1 - xMin) * WM2min)),
                sqrt(min(Q2max, xMax * yMax * sM2, 1000 * WM2max))]


def _read_compass23_table(path, start, n_rows):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    for line in lines[start:start + n_rows]:
        line = line.replace("[", "").replace("]", "").replace(";", ",").replace(",", " ")
        rows.append([float(v) for v in line.split()])
    return rows


def _add_compass23_points(ds, rows, diff_var, h_2):
    for i, row in enumerate(rows):
        if diff_var == "x":
            x_bin, z_bin, pT_bin = (row[0], row[1]), (0.2, 1.), (0.1, 1.6)
        elif diff_var == "z":
            x_bin, z_bin, pT_bin = (0.003, 0.7), (row[0], row[1]), (0.1, 1.6)
        else:  # "pt"
            x_bin, z_bin, pT_bin = (0.003, 0.7), (0.2, 1.), (row[0], row[1])

        Q_min, Q_max = Qbounds(x_bin[0], x_bin[1])

        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=12001,
            ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=h_2, proc_id_weight=proc_id_weight,
            s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[7]),
            x_min=x_bin[0], x_max=x_bin[1], x_avg=row[2],
            z_min=z_bin[0], z_max=z_bin[1], z_avg=row[4],
            pT_min=pT_bin[0], pT_max=pT_bin[1], pT_avg=row[5],
            M_target=M_proton, M_product=m_pion,
            includeCuts=includeCuts,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            thFactor=1.,
            xSec=row[8],
            uncorrErr_0=row[10],
        ))


#%% -- h+, dx
rows = _read_compass23_table(path_to_data + "Sivers.dat", 6, 9)
ds = DataSet.empty("SIDIS", name="compass23.sivers.h+.dx", comment="COMPASS SSA-Sivers h+ (differential in x)",
                    reference="2401.00309", normErr=[], isNormalized=False)
_add_compass23_points(ds, rows, "x", 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, dz
rows = _read_compass23_table(path_to_data + "Sivers.dat", 19, 8)
ds = DataSet.empty("SIDIS", name="compass23.sivers.h+.dz", comment="COMPASS SSA-Sivers h+ (differential in z)",
                    reference="2401.00309", normErr=[], isNormalized=False)
_add_compass23_points(ds, rows, "z", 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, dpt
rows = _read_compass23_table(path_to_data + "Sivers.dat", 31, 9)
ds = DataSet.empty("SIDIS", name="compass23.sivers.h+.dpt", comment="COMPASS SSA-Sivers h+ (differential in pt)",
                    reference="2401.00309", normErr=[], isNormalized=False)
_add_compass23_points(ds, rows, "pt", 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, dx
rows = _read_compass23_table(path_to_data + "Sivers.dat", 44, 9)
ds = DataSet.empty("SIDIS", name="compass23.sivers.h-.dx", comment="COMPASS SSA-Sivers h- (differential in x)",
                    reference="2401.00309", normErr=[], isNormalized=False)
_add_compass23_points(ds, rows, "x", -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, dz
rows = _read_compass23_table(path_to_data + "Sivers.dat", 57, 8)
ds = DataSet.empty("SIDIS", name="compass23.sivers.h-.dz", comment="COMPASS SSA-Sivers h- (differential in z)",
                    reference="2401.00309", normErr=[], isNormalized=False)
_add_compass23_points(ds, rows, "z", -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, dpt
rows = _read_compass23_table(path_to_data + "Sivers.dat", 69, 9)
ds = DataSet.empty("SIDIS", name="compass23.sivers.h-.dpt", comment="COMPASS SSA-Sivers h- (differential in pt)",
                    reference="2401.00309", normErr=[], isNormalized=False)
_add_compass23_points(ds, rows, "pt", -12)
ds.save_csv(path_to_save + ds.name + ".csv")
