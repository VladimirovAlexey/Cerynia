"""
Parsing of COMPASS (2008 publication) transverse single-spin asymmetry
A_UT^Sivers for pi+/pi-/k+/k-/k0, deuteron target: "all hadrons" (main)
and "leading hadrons" samples.

Source: /data/arTeMiDe_Repository/data/COMPASS/0802.2160/durham_0802.2160.txt
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from the "compass08.sivers.*" (Sivers) blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSSA_Compass08.py -- s, x/z/pT/Q
bins, cuts, and error treatment kept identical to the old parsing. The
sibling Collins-asymmetry datasets (compass08.collins.*) built by the same
old script are NOT ported -- Sivers only, per user instruction.

NOTE: process code. Old parsing used [1,1,h_2,12002] (primary, h_2: pi+=1,
pi-=-1, k+=2, k-=-2, k0=3) / [1,1,h_2,2002] (weight), with h_1=1 (proton)
despite the raw file's own title "Collins and Sivers asymmetries for pions
and kaons in muon-deuteron DIS". Per user confirmation, updated to primary
[1,12,h_2,12001] / weight [1,12,h_2,2001] -- h_1=12 (deuteron); h_2: pi+=1,
pi-=-1, k+=2, k-=-2, k0=20 (NOTE: k0=20 here, NOT the same as k0=3 used in
the old scheme or the h_2=12 isoscalar-hadron code used for COMPASS
h+/h- -- k0 gets its own new code). proc_id=12001/weight=2001 match the
now-uniform SIDIS-Sivers codes seen in JLab.py/hermes09.py/hermes3D.py/
COMPASS23.py.

NOTE: "leading hadron" datasets. The old script also builds 15 more
datasets ("...leading" blocks, e.g. compass08.sivers.pi+leading.dx) with
active SaveToCSV calls, but NONE of these exist in the actual old on-disk
DataLib/Sivers/ (only the 15 non-leading ones do) -- and their old
proc_current=[1,1,12001] was only 3 elements (missing h_2 entirely),
looking unfinished. Per user instruction, these ARE ported here (with the
same corrected process code as the main sets, h_1=12/proc_id=12001/
weight=2001), but as a fully separate physical sample (different z_min
cut: 0.25 instead of 0.2, per old parsing) -- NOT a duplicate/subset of
the main sets.

NOTE: naming convention, given by the user ("place .leading at the end"):
old "compass08.sivers.pi+leading.dx" -> "compass08.sivers.pi+.dx.leading"
(moving ".leading" to the very end of the name, after the differential
variable, rather than gluing it to the hadron label as the old scheme did).

NOTE: raw table columns reorder depending on which variable is
differential: the first 3 fields are always [avg, low, high] of the
binned variable; Q2 and y follow; then the OTHER two kinematic variables'
averages (order swaps between dx/dz/dpt blocks, see _add_compass08_points);
then Collins value/error (unused, Sivers-only per user instruction), then
Sivers value/error (used as xSec/uncorrErr_0).

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention).

Datasets, 30 total (arXiv:0802.2160):
    compass08.sivers.{pi+,pi-,k+,k-,k0}.{dx,dz,dpt}          -- 15 main sets
    compass08.sivers.{pi+,pi-,k+,k-,k0}.{dx,dz,dpt}.leading  -- 15 leading sets
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/COMPASS/0802.2160/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_proton = 0.938
m_pion = 0.139
m_kaon = 0.494

### 160 GeV beam on a fixed deuteron target (as in old parsing)
s = 2. * 160. * 0.938 + 0.938**2
includeCuts = True
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.9, 25., 10000.]

### deuteron target (per user instruction)
ps_def, h_1 = 1, 12
proc_id, proc_id_weight = 12001, 2001

### fixed bin-integration ranges (as in old parsing); z_min differs main/leading
X_FIXED = (0.03, 1.)
PT_FIXED = (0.1, 10.)
Z_FIXED_MAIN = (0.2, 1.)
Z_FIXED_LEADING = (0.25, 1.)


def Qbounds(xMin, xMax):
    """
    Kinematic Q bounds given an x window (ported verbatim from old parsing,
    same as COMPASS23.py):
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


def _read_compass08_table(path, start, n_rows):
    with open(path) as f:
        lines = [line.rstrip() for line in f]
    ### old parsing does `del data_from_f[0:13]` once up front (strips the file
    ### header/title block); every block offset below is relative to that
    lines = lines[13:]
    rows = []
    for line in lines[start:start + n_rows]:
        line = line.replace("[", "").replace("]", "").replace(",", " ")
        rows.append([float(v) for v in line.split()])
    return rows


def _add_compass08_points(ds, rows, diff_var, h_2, z_fixed, m_product):
    """
    diff_var: which variable is binned per-row ("x", "z", or "pt"); row
    layout is [avg, low, high, Q2, y(unused), other1_avg, other2_avg,
    Col_val(unused), Col_err(unused), Siv_val, Siv_err], where
    (other1, other2) = (z, pT) for diff_var="x", (x, pT) for "z", (x, z)
    for "pt" -- matching old parsing's column order exactly.
    """
    for i, row in enumerate(rows):
        if diff_var == "x":
            x_bin, z_bin, pT_bin = (row[1], row[2]), z_fixed, PT_FIXED
        elif diff_var == "z":
            x_bin, z_bin, pT_bin = X_FIXED, (row[1], row[2]), PT_FIXED
        else:  # "pt"
            x_bin, z_bin, pT_bin = X_FIXED, z_fixed, (row[1], row[2])

        Q_min, Q_max = Qbounds(x_bin[0], x_bin[1])

        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=h_2, proc_id_weight=proc_id_weight,
            s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[3]),
            x_min=x_bin[0], x_max=x_bin[1],
            z_min=z_bin[0], z_max=z_bin[1],
            pT_min=pT_bin[0], pT_max=pT_bin[1],
            M_target=M_proton, M_product=m_product,
            includeCuts=includeCuts,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            thFactor=1.,
            xSec=row[9],
            uncorrErr_0=row[10],
        ))


#%% -- pi+, dx
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 0, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi+.dx", comment="COMPASS SSA-Sivers pi+ (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", 1, Z_FIXED_MAIN, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, dpt
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 11, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi+.dpt", comment="COMPASS SSA-Sivers pi+ (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", 1, Z_FIXED_MAIN, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, dz
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 22, 8)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi+.dz", comment="COMPASS SSA-Sivers pi+ (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", 1, Z_FIXED_MAIN, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, dx
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 34, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi-.dx", comment="COMPASS SSA-Sivers pi- (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", -1, Z_FIXED_MAIN, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, dpt
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 45, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi-.dpt", comment="COMPASS SSA-Sivers pi- (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", -1, Z_FIXED_MAIN, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, dz
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 56, 8)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi-.dz", comment="COMPASS SSA-Sivers pi- (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", -1, Z_FIXED_MAIN, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, dx
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 68, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k+.dx", comment="COMPASS SSA-Sivers k+ (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", 2, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, dpt
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 79, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k+.dpt", comment="COMPASS SSA-Sivers k+ (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", 2, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, dz
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 90, 8)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k+.dz", comment="COMPASS SSA-Sivers k+ (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", 2, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, dx
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 102, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k-.dx", comment="COMPASS SSA-Sivers k- (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", -2, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, dpt
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 113, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k-.dpt", comment="COMPASS SSA-Sivers k- (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", -2, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, dz
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 124, 8)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k-.dz", comment="COMPASS SSA-Sivers k- (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", -2, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k0, dx
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 136, 5)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k0.dx", comment="COMPASS SSA-Sivers k0 (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", 20, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k0, dpt
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 143, 5)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k0.dpt", comment="COMPASS SSA-Sivers k0 (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", 20, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k0, dz
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 150, 6)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k0.dz", comment="COMPASS SSA-Sivers k0 (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", 20, Z_FIXED_MAIN, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, dx, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 162, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi+.dx.leading", comment="COMPASS SSA-Sivers pi+ leading (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", 1, Z_FIXED_LEADING, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, dpt, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 173, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi+.dpt.leading", comment="COMPASS SSA-Sivers pi+ leading (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", 1, Z_FIXED_LEADING, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, dz, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 184, 7)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi+.dz.leading", comment="COMPASS SSA-Sivers pi+ leading (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", 1, Z_FIXED_LEADING, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, dx, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 195, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi-.dx.leading", comment="COMPASS SSA-Sivers pi- leading (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", -1, Z_FIXED_LEADING, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, dpt, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 206, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi-.dpt.leading", comment="COMPASS SSA-Sivers pi- leading (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", -1, Z_FIXED_LEADING, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, dz, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 217, 7)
ds = DataSet.empty("SIDIS", name="compass08.sivers.pi-.dz.leading", comment="COMPASS SSA-Sivers pi- leading (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", -1, Z_FIXED_LEADING, m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, dx, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 228, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k+.dx.leading", comment="COMPASS SSA-Sivers k+ leading (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", 2, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, dpt, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 239, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k+.dpt.leading", comment="COMPASS SSA-Sivers k+ leading (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", 2, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, dz, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 250, 7)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k+.dz.leading", comment="COMPASS SSA-Sivers k+ leading (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", 2, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, dx, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 261, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k-.dx.leading", comment="COMPASS SSA-Sivers k- leading (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", -2, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, dpt, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 272, 9)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k-.dpt.leading", comment="COMPASS SSA-Sivers k- leading (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", -2, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, dz, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 283, 7)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k-.dz.leading", comment="COMPASS SSA-Sivers k- leading (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", -2, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k0, dx, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 294, 5)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k0.dx.leading", comment="COMPASS SSA-Sivers k0 leading (differential in x)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "x", 20, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k0, dpt, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 301, 5)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k0.dpt.leading", comment="COMPASS SSA-Sivers k0 leading (differential in pt)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "pt", 20, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k0, dz, leading
rows = _read_compass08_table(path_to_data + "durham_0802.2160.txt", 308, 5)
ds = DataSet.empty("SIDIS", name="compass08.sivers.k0.dz.leading", comment="COMPASS SSA-Sivers k0 leading (differential in z)",
                    reference="0802.2160", normErr=[], isNormalized=False)
_add_compass08_points(ds, rows, "z", 20, Z_FIXED_LEADING, m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")
