"""
Parsing of COMPASS (2016 publication, 2010 proton-target data) double-spin
asymmetry A_LT^cos(phi_h-phi_S) for isoscalar h+/h-, proton target, 3
z-ranges x 3 differential variables, joined across 4 Q-windows.

Source: /data/arTeMiDe_Repository/data/COMPASS/1609.07374/ALT_cosHmS/Zgt*_{hp,hm}_{pt,x,z}_ALT_cosHmS.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/wgt/*.csv

Ported from DataProcessor/OtherPrograms/DataParsing/ReadALT_Compass.py --
s, x/z/pT/Q bins, cuts, and error treatment kept identical to the old
parsing. Structurally the same measurement/binning scheme as
Sivers/COMPASS16.py (same paper, same z-ranges/Q-windows/differential
variables), but a DIFFERENT observable (A_LT, not A_Sivers) and a
DIFFERENT target.

NOTE: target is PROTON, not deuteron. The raw files' own header literally
says "COMPASS proton 2010 data" for every one of the 18 files -- confirmed
with the user before writing this (initially assumed deuteron, by analogy
with Sivers/COMPASS16.py, since both share the same 1609.07374 reference).
User confirmed: this measurement really is proton-target (the "2010"
in the raw header is the data-collection year; "compass16" in the dataset
name refers to the 2016 publication year, not a distinct deuteron run).
So h_1=1 here, unlike every other COMPASS Sivers/wgt sub-case in this
migration (COMPASS08/16/23 Sivers are all deuteron, h_1=12).

NOTE: process code, confirmed by the user. Primary [1,1,h_2,13001], weight
[1,1,h_2,2001] (h_2: 12=h+, -12=h- -- same isoscalar-hadron code as the
deuteron-target COMPASS cases; confirmed this code is purely a hadron-type
tag, not target-dependent). proc_id=13001 is new here (distinguishing the
A_LT/wgt observable from Sivers' 12001); weight proc_id=2001 is the SAME
weighting/normalization process as every other SIDIS Sivers/wgt sub-case
in this migration -- the underlying unpolarized cross section used to
normalize an asymmetry doesn't depend on which spin observable is in the
numerator.

NOTE: old library category/dataset-name segment "ALT" renamed to "wgt"
(matching the Cerynia category folder name), per user instruction:
"compass16.ALT.h+.1<z<2.dpt" -> "compass16.wgt.h+.1<z<2.dpt.dQ".

NOTE: "joined across Q-windows" (.dQ suffix), per user instruction --
same treatment as Sivers/COMPASS16.py. Unlike the Sivers case, THIS old
script already saves both the joined ("compass16.ALT.<...>", no Q-window
in the name) AND per-Q-window ("...binQ1" through "binQ4") datasets
directly (no need to infer the joined form from a second script -- one
script builds all 5 variants explicitly, appending each point to both the
joined DataSet and its matching per-Q-window DataSet). Only the already-
joined form is ported here, renamed with ".dQ"; the 4 binQ1-4 variants are
pure subsets, not ported separately (same pattern as Sivers/COMPASS16.py
and E615's dropped per-bin child datasets).

NOTE: raw table column layout matches Sivers/COMPASS16.py exactly: columns
0,1 = bin edges of the differential variable; column 2 = x_avg, column 4 =
z_avg, column 5 = pT_avg (fixed positions regardless of which variable is
differential); column 7 = Q2_avg; column 8 = xSec (Asym); column 10 =
uncorrErr (er_tot; column 9, er_stat, is skipped). For whichever variable
is non-differential and non-Q: x gets xbounds(Q2min,Q2max) (ported
verbatim, same formula as Sivers/COMPASS16.py); z and pT get fixed
literals (z: per z-range group; pT: always [0.1,10.]).

NOTE: normErr=[0.03], a dilution factor ("delution factor according to
Bakur" per old code's comment) -- ported verbatim, a real physics value
present in the old script, not something introduced here.

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention).

Datasets, 18 total (arXiv:1609.07374), each double-differential in Q (4
windows: 1-2, 2-2.5, 2.5-4, 4-9 GeV) and the named variable:
    compass16.wgt.{h+,h-}.{1<z<2,1<z,2<z}.{dpt,dx,dz}
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/COMPASS/1609.07374/ALT_cosHmS/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/wgt/"

M_proton = 0.938
m_pion = 0.139

### 160 GeV beam on a fixed proton target (as in old parsing)
s = 2. * 160. * 0.938 + 0.938**2
includeCuts = True
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.9, 10., 10000.]
### dilution factor (as in old parsing, "according to Bakur")
normErr = [0.03]

### proton target (confirmed by the user, unlike the deuteron-target COMPASS Sivers cases)
ps_def, h_1 = 1, 1
proc_id, proc_id_weight = 13001, 2001

PT_FIXED = (0.1, 10.)


def xbounds(Q2min, Q2max):
    """
    Kinematic x bounds given a Q^2 window (ported verbatim from old parsing,
    same as Sivers/COMPASS16.py):
    xmin = MAX{xmin, Q2/(yMax*(s-M^2)), Q2/(Q2+W2max-M^2)}
    xmax = MIN{xmax, Q2/(yMin*(s-M^2)), Q2/(Q2+W2min-M^2)}
    """
    xmin, xmax = 0.003, 0.9
    WM2min, WM2max = 10. - 0.938**2, 10000. - 0.938**2
    yMin, yMax = 0.1, 0.9
    sM2 = 2. * 160. * 0.938
    return [max(xmin, Q2min / (yMax * sM2), Q2min / (Q2min + WM2max)),
            min(xmax, Q2max / (yMin * sM2), Q2max / (Q2max + WM2min))]


def _read_compass16_alt_table(path):
    """
    Read a COMPASS16 ALT table: blocks separated by "Q^{2}/(GeV/c)^2"
    marker lines, data rows starting with "[low;high]\\t...". Returns a
    list of (Q2_min, Q2_max, rows).
    """
    with open(path) as f:
        lines = [line.rstrip() for line in f]
    blocks = []
    Q2_current = None
    rows = None
    for line in lines[5:]:
        if line == "" or "<x>" in line:
            continue
        if "Q^{2}/(GeV/c)^{2}" in line:
            if rows is not None:
                blocks.append((Q2_current[0], Q2_current[1], rows))
            if "1<Q^{2}/(GeV/c)^{2}<4" in line:
                Q2_current = (1., 4.)
            elif "4<Q^{2}/(GeV/c)^{2}<6.25" in line:
                Q2_current = (4., 6.25)
            elif "6.25<Q^{2}/(GeV/c)^{2}<16" in line:
                Q2_current = (6.25, 16.)
            elif "16<Q^{2}/(GeV/c)^{2}<81" in line:
                Q2_current = (16., 81.)
            else:
                raise ValueError(f"CANNOT DETERMINE Q in line: {line}")
            rows = []
        elif line[0] == "[":
            parts = line.replace("[", "").replace("]", "").replace(";", " ").split()
            rows.append([float(v) for v in parts])
        else:
            raise ValueError(f"DO NOT UNDERSTAND the line: {line}")
    if rows is not None:
        blocks.append((Q2_current[0], Q2_current[1], rows))
    return blocks


def _add_compass16_alt_points(ds, path, var, z_bin, h_2):
    """var: "pt", "x", or "z" -- which variable is differential (bin edges from row[0],row[1])."""
    idx = 0
    for Q2_min, Q2_max, rows in _read_compass16_alt_table(path):
        Q_min, Q_max = sqrt(Q2_min), sqrt(Q2_max)
        x_kin = xbounds(Q2_min, Q2_max)
        for row in rows:
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
ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.1<z<2.dpt", comment="COMPASS16 DSSA-A_LT h- from p, 0.1<z<0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01ls02_hm_pt_ALT_cosHmS.dat", "pt", (0.1, 0.2), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.1<z<2.dx", comment="COMPASS16 DSSA-A_LT h- from p, 0.1<z<0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01ls02_hm_x_ALT_cosHmS.dat", "x", (0.1, 0.2), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.1<z<2.dz", comment="COMPASS16 DSSA-A_LT h- from p, 0.1<z<0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01ls02_hm_z_ALT_cosHmS.dat", "z", (0.1, 0.2), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, 1<z (z>0.1)
ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.1<z.dpt", comment="COMPASS16 DSSA-A_LT h- from p, z>0.1 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01_hm_pt_ALT_cosHmS.dat", "pt", (0.1, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.1<z.dx", comment="COMPASS16 DSSA-A_LT h- from p, z>0.1 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01_hm_x_ALT_cosHmS.dat", "x", (0.1, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.1<z.dz", comment="COMPASS16 DSSA-A_LT h- from p, z>0.1 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01_hm_z_ALT_cosHmS.dat", "z", (0.1, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h-, 2<z (z>0.2)
ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.2<z.dpt", comment="COMPASS16 DSSA-A_LT h- from p, z>0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt02_hm_pt_ALT_cosHmS.dat", "pt", (0.2, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.2<z.dx", comment="COMPASS16 DSSA-A_LT h- from p, z>0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt02_hm_x_ALT_cosHmS.dat", "x", (0.2, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h-.2<z.dz", comment="COMPASS16 DSSA-A_LT h- from p, z>0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt02_hm_z_ALT_cosHmS.dat", "z", (0.2, 1.), -12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, 1<z<2 (0.1<z<0.2)
ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.1<z<2.dpt", comment="COMPASS16 DSSA-A_LT h+ from p, 0.1<z<0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01ls02_hp_pt_ALT_cosHmS.dat", "pt", (0.1, 0.2), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.1<z<2.dx", comment="COMPASS16 DSSA-A_LT h+ from p, 0.1<z<0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01ls02_hp_x_ALT_cosHmS.dat", "x", (0.1, 0.2), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.1<z<2.dz", comment="COMPASS16 DSSA-A_LT h+ from p, 0.1<z<0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01ls02_hp_z_ALT_cosHmS.dat", "z", (0.1, 0.2), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, 1<z (z>0.1)
ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.1<z.dpt", comment="COMPASS16 DSSA-A_LT h+ from p, z>0.1 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01_hp_pt_ALT_cosHmS.dat", "pt", (0.1, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.1<z.dx", comment="COMPASS16 DSSA-A_LT h+ from p, z>0.1 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01_hp_x_ALT_cosHmS.dat", "x", (0.1, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.1<z.dz", comment="COMPASS16 DSSA-A_LT h+ from p, z>0.1 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt01_hp_z_ALT_cosHmS.dat", "z", (0.1, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- h+, 2<z (z>0.2)
ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.2<z.dpt", comment="COMPASS16 DSSA-A_LT h+ from p, z>0.2 (differential in pt, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt02_hp_pt_ALT_cosHmS.dat", "pt", (0.2, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.2<z.dx", comment="COMPASS16 DSSA-A_LT h+ from p, z>0.2 (differential in x, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt02_hp_x_ALT_cosHmS.dat", "x", (0.2, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")

ds = DataSet.empty("SIDIS", name="compass16.wgt.h+.2<z.dz", comment="COMPASS16 DSSA-A_LT h+ from p, z>0.2 (differential in z, joined across 4 Q-windows)",
                    reference="1609.07374", normErr=normErr, isNormalized=False)
_add_compass16_alt_points(ds, path_to_data + "Zgt02_hp_z_ALT_cosHmS.dat", "z", (0.2, 1.), 12)
ds.save_csv(path_to_save + ds.name + ".csv")
