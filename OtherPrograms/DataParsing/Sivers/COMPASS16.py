"""
Parsing of COMPASS (2016 publication, 2010 proton-target data) transverse
single-spin asymmetry A_UT^Sivers for isoscalar h+/h-, proton target, 3
z-ranges x 3 differential variables, joined across 4 Q-windows.

Source: /data/arTeMiDe_Repository/data/COMPASS/1609.07374/Zgt*_{n,p}_{pt,x,z}_Siv.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from ReadSSA_Compass16_v2.py's "full" blocks (NOT the older
ReadSSA_Compass16.py, nor the "binQ1"-"binQ4" split from
ReadSSA_Compass16_v2_additional.py -- see NOTE below on the investigation
that resolved this).

NOTE: script investigation, done before writing this file. Three old
scripts exist for COMPASS16:
  - ReadSSA_Compass16.py (oldest): produces 72 datasets named e.g.
    "compass16.sivers.h-.1<z<2.1<Q<2.dpt", with the Q-window spelled out in
    the name. Compared byte-for-byte against the newer scripts' equivalent
    block: identical raw file, identical line slice, identical kinematics
    -- but its proc_current was only 3 elements (missing h_2 entirely, a
    completeness bug). Confirmed superseded; not used at all.
  - ReadSSA_Compass16_v2_additional.py: produces 72 datasets named
    "compass16.sivers.h-.1<z<2.dpt.binQ1" (through binQ4) -- same 4
    Q-windows as above, renamed, with a complete proc_current. Each binQ
    file holds the SAME points that appear (as a subset) inside v2.py's
    "full" dataset for that (charge, z-range, variable) -- see next point.
  - ReadSSA_Compass16_v2.py: produces 18 "full" datasets (e.g.
    "compass16.sivers.h-.1<z<2.dpt", no Q-window in the name) by looping
    over the SAME 4 Q-windows internally and appending all their points
    into ONE DataSet object, calling SaveToCSV only once at the end (the
    intermediate SaveToCSV calls for each Q-sub-block are commented out).
    I.e. v2.py's "full" sets are ALREADY the union of the 4 binQ1-4 sets --
    not a genuinely Q-integrated single measurement.
  Per user instruction ("join these separate binQ into set with .dQ
  suffix"), this file ports v2.py's already-joined 18 datasets directly
  (each spanning all 4 Q-windows, double-differential in Q and pt/x/z) and
  renames them with a ".dQ" suffix to make the multi-Q-window nature
  explicit. The 72 binQ1-4 datasets are NOT ported separately -- they are
  pure subsets of these 18 (same "redundant per-bin child dataset" pattern
  as E615's dropped dQ-4.05.../dxF-0.0... children).

NOTE: process code. Old parsing used [1,1,h_2,12101] (primary) /
[1,1,h_2,2101] (weight). Initially "corrected" h_1 to 12 (deuteron) here,
by analogy with the deuteron-target COMPASS08/23 Sivers cases -- WRONG.
This measurement (and the sibling wgt/COMPASS16.py A_LT data, same raw
data folder /COMPASS/1609.07374/) is confirmed PROTON target: the paper
itself (Phys.Lett.B770(2017)138, arXiv:1609.07374; PDF in the raw data
folder) has the abstract "Eight PROTON transverse-spin-dependent azimuthal
asymmetries are extracted... from the COMPASS 2010 semi-inclusive hadron
measurements... for the azimuthal asymmetries induced by the Sivers
[TMD]..." -- i.e. this Sivers-2D dataset and the wgt A_LT dataset are the
SAME underlying COMPASS 2010 proton-target measurement, just different
asymmetries extracted from it, not two different target/runs. Confirmed
and fixed per user instruction: h_1=1 (proton), in both the primary AND
weight process. Final: primary [1,1,h_2,12001] / weight [1,1,h_2,2001]
(h_2: 12=h+, -12=h- -- this hadron-type code is target-independent,
confirmed separately).

NOTE: raw table column layout is fixed regardless of which variable is
differential: columns 0,1 are always the bin edges of whichever variable
is differential; column 2 = x_avg, column 4 = z_avg, column 5 = pT_avg,
column 7 = Q2_avg, column 8 = xSec (Sivers value), column 10 = uncorrErr
(Sivers total error) -- column 9 (Sivers stat error) is skipped, same as
column-skipping pattern in COMPASS08/23. For whichever variable is NOT
differential and NOT the binned Q-window: x gets a kinematically-derived
bound via xbounds(Q2min,Q2max) (ported verbatim, x is tightly correlated
with Q so a fixed literal wouldn't make sense per Q-window); z and pT get
simple FIXED literal ranges (z: per z-range group, see Z_GROUPS; pT:
always [0.1,10.]).

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention). M_product
uses the pion mass as a stand-in for the unidentified isoscalar hadron
(same convention as compass.d.h+/h- and COMPASS08/23).

Datasets, 18 total (arXiv:1609.07374), each double-differential in Q (4
windows: 1-2, 2-2.5, 2.5-4, 4-9 GeV) and the named variable:
    compass16.sivers.{h+,h-}.{1<z<2,1<z,2<z}.{dpt,dx,dz}.dQ
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/COMPASS/1609.07374/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_proton = 0.938
m_pion = 0.139

### 160 GeV beam on a fixed proton target (as in old parsing)
s = 2. * 160. * 0.938 + 0.938**2
includeCuts = True
### y, W^2 cuts (as in old parsing; W2min=10, NOT 25 like COMPASS08/23)
cutParams = [0.1, 0.9, 10., 10000.]

### proton target (confirmed by the user; see module NOTE -- this measurement
### is the same COMPASS 2010 proton run as wgt/COMPASS16.py, not deuteron)
ps_def, h_1 = 1, 1
proc_id, proc_id_weight = 12001, 2001

PT_FIXED = (0.1, 10.)
Q2_WINDOWS = [(1., 4.), (4., 6.25), (6.25, 16.), (16., 81.)]
### (start, n_rows) per Q-window, 0-indexed into the raw file's line list
PT_SLICES = [(6, 5), (15, 5), (24, 5), (33, 5)]
X_SLICES = [(6, 7), (17, 6), (27, 7), (38, 6)]
Z_SLICES = [(6, 5), (15, 5), (24, 5), (33, 5)]


def xbounds(Q2min, Q2max):
    """
    Kinematic x bounds given a Q^2 window (ported verbatim from old parsing):
    xmin = MAX{xmin, Q2/(yMax*(s-M^2)), Q2/(Q2+W2max-M^2)}
    xmax = MIN{xmax, Q2/(yMin*(s-M^2)), Q2/(Q2+W2min-M^2)}
    """
    xmin, xmax = 0.003, 0.9
    WM2min, WM2max = 10. - 0.938**2, 10000. - 0.938**2
    yMin, yMax = 0.1, 0.9
    sM2 = 2. * 160. * 0.938
    return [max(xmin, Q2min / (yMax * sM2), Q2min / (Q2min + WM2max)),
            min(xmax, Q2max / (yMin * sM2), Q2max / (Q2max + WM2min))]


def _read_compass16_rows(path, start, n_rows):
    with open(path) as f:
        lines = [line.rstrip() for line in f]
    rows = []
    for line in lines[start:start + n_rows]:
        line = line.replace("[", "").replace("]", "").replace(";", " ")
        rows.append([float(v) for v in line.split()])
    return rows


def _add_compass16_points(ds, path, var, z_bin, h_2):
    """
    var: "pt", "x", or "z" -- which variable is differential (bin edges
    from columns 0,1). Joins all 4 Q-windows into ds.
    """
    slices = {"pt": PT_SLICES, "x": X_SLICES, "z": Z_SLICES}[var]
    idx = 0
    for (start, n_rows), (Q2_min, Q2_max) in zip(slices, Q2_WINDOWS):
        Q_min, Q_max = sqrt(Q2_min), sqrt(Q2_max)
        x_kin = xbounds(Q2_min, Q2_max)
        for row in _read_compass16_rows(path, start, n_rows):
            if var == "pt":
                pT_bin, x_bin, z_b = (row[0], row[1]), x_kin, z_bin
            elif var == "x":
                pT_bin, x_bin, z_b = PT_FIXED, (row[0], row[1]), z_bin
            else:  # "z"
                pT_bin, x_bin, z_b = PT_FIXED, x_kin, (row[0], row[1])

            ds.add_point(dict(
                id=f"{ds.name}.{idx}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
                ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=h_2, proc_id_weight=proc_id_weight,
                s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[7]),
                x_min=x_bin[0], x_max=x_bin[1], x_avg=row[2],
                z_min=z_b[0], z_max=z_b[1], z_avg=row[4],
                pT_min=pT_bin[0], pT_max=pT_bin[1], pT_avg=row[5],
                M_target=M_proton, M_product=m_pion,
                includeCuts=includeCuts,
                cutParams_0=cutParams[0], cutParams_1=cutParams[1],
                cutParams_2=cutParams[2], cutParams_3=cutParams[3],
                thFactor=1.,
                xSec=row[8],
                uncorrErr_0=row[10],
            ))
            idx += 1


#%% -- h-, 1<z<2 (0.1<z<0.2)
ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.1<z<2.dpt.dQ", comment="COMPASS16 SSA-Sivers h-, 0.1<z<0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01ls02_n_pt_Siv.dat", "pt", (0.1, 0.2), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.1<z<2.dx.dQ", comment="COMPASS16 SSA-Sivers h-, 0.1<z<0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01ls02_n_x_Siv.dat", "x", (0.1, 0.2), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.1<z<2.dz.dQ", comment="COMPASS16 SSA-Sivers h-, 0.1<z<0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01ls02_n_z_Siv.dat", "z", (0.1, 0.2), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, 1<z (z>0.1)
ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.1<z.dpt.dQ", comment="COMPASS16 SSA-Sivers h-, z>0.1 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01_n_pt_Siv.dat", "pt", (0.1, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.1<z.dx.dQ", comment="COMPASS16 SSA-Sivers h-, z>0.1 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01_n_x_Siv.dat", "x", (0.1, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.1<z.dz.dQ", comment="COMPASS16 SSA-Sivers h-, z>0.1 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01_n_z_Siv.dat", "z", (0.1, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, 2<z (z>0.2)
ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.2<z.dpt.dQ", comment="COMPASS16 SSA-Sivers h-, z>0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt02_n_pt_Siv.dat", "pt", (0.2, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.2<z.dx.dQ", comment="COMPASS16 SSA-Sivers h-, z>0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt02_n_x_Siv.dat", "x", (0.2, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h-.2<z.dz.dQ", comment="COMPASS16 SSA-Sivers h-, z>0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt02_n_z_Siv.dat", "z", (0.2, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, 1<z<2 (0.1<z<0.2)
ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.1<z<2.dpt.dQ", comment="COMPASS16 SSA-Sivers h+, 0.1<z<0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01ls02_p_pt_Siv.dat", "pt", (0.1, 0.2), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.1<z<2.dx.dQ", comment="COMPASS16 SSA-Sivers h+, 0.1<z<0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01ls02_p_x_Siv.dat", "x", (0.1, 0.2), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.1<z<2.dz.dQ", comment="COMPASS16 SSA-Sivers h+, 0.1<z<0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01ls02_p_z_Siv.dat", "z", (0.1, 0.2), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, 1<z (z>0.1)
ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.1<z.dpt.dQ", comment="COMPASS16 SSA-Sivers h+, z>0.1 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01_p_pt_Siv.dat", "pt", (0.1, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.1<z.dx.dQ", comment="COMPASS16 SSA-Sivers h+, z>0.1 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01_p_x_Siv.dat", "x", (0.1, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.1<z.dz.dQ", comment="COMPASS16 SSA-Sivers h+, z>0.1 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt01_p_z_Siv.dat", "z", (0.1, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, 2<z (z>0.2)
ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.2<z.dpt.dQ", comment="COMPASS16 SSA-Sivers h+, z>0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt02_p_pt_Siv.dat", "pt", (0.2, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.2<z.dx.dQ", comment="COMPASS16 SSA-Sivers h+, z>0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt02_p_x_Siv.dat", "x", (0.2, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.sivers.h+.2<z.dz.dQ", comment="COMPASS16 SSA-Sivers h+, z>0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=[], isNormalized=False)
_add_compass16_points(ds, path_to_data + "Zgt02_p_z_Siv.dat", "z", (0.2, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")
