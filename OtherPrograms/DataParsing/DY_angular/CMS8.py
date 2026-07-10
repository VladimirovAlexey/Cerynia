"""
Parsing of CMS 8 TeV Z-boson angular coefficients A0-A4 (Collins-Soper frame),
2 rapidity bins, near the Z peak.

Source: /data/arTeMiDe_Repository/data/CMS/1504.03512/y0to1.csv, y1to21.csv
        (HEPData record for arXiv:1504.03512)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY_angular/*.csv

Ported from the "CMS8-A*" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(angularCMS).py --
definitions (process code, s, Q-range, y-ranges, cuts, error split) are kept
identical to the old parsing. See ATLAS8.py (same DY_angular category) for
the shared schema/weight-column notes.

NOTE: proc_id numbering (20-27 for A0-A7) is the same convention as ATLAS8.py,
confirmed correct by the user.

NOTE: the old library's "CMS8-Auu-*" datasets are NOT ported -- same bug
pattern found in ATLAS8.py's dropped "A8_Auu_*" sets. Old code built
CMS8-Auu-00y10/10y21 by reading the *same* raw-file line range as
CMS8-A0-00y10/10y21 (both `data_from_f[15:23]`), i.e. Auu was byte-for-byte
identical to A0, just tagged proc_id=3 instead of 20. Per user instruction
("weighting processes is the same" as ATLAS8), every point here instead
carries the weight-process columns (ps_def_weight, h_1_weight, h_2_weight,
proc_id_weight) = (1,1,1,3), same mechanism as ATLAS8.py.

NOTE: thFactor. Old parsing used thFactor=2/(qT_max-qT_min) here (a bin-width
-dividing, y-symmetrizing factor, as if this were an ordinary differential
cross section) -- inconsistent with ATLAS8.py's angular thFactor=1 for the
same kind of ratio observable. Per user decision, set thFactor=1 here too,
treating the old CMS8 value as a leftover from copy-pasting cross-section
parsing code rather than an intentional difference.

Datasets (arXiv:1504.03512), 10 total (A0-A4 x 2 y-bins; no A5/A6/A7 measured
in this paper):
    CMS8-A0-00y10, CMS8-A1-00y10, CMS8-A2-00y10, CMS8-A3-00y10, CMS8-A4-00y10
    CMS8-A0-10y21, CMS8-A1-10y21, CMS8-A2-10y21, CMS8-A3-10y21, CMS8-A4-10y21
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/CMS/1504.03512/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY_angular/"

M_Z = 91.1876  # forces <Q> = Z-boson mass, matching old parsing's explicit p["<Q>"]=M_Z

### 8 TeV, near the Z peak (as in old parsing)
s = 8000.**2
Q_min, Q_max = 81., 101.
### CMS 8 lepton fiducial cuts, unused since includeCuts=False (as in old parsing)
includeCuts = False
cutParams = [20., 20., 2., 4.5]
### weighting process: plain unpolarized DY cross section (ps_def,h_1,h_2,proc_id),
### same mechanism as ATLAS8.py, per user confirmation
ps_def_weight, h_1_weight, h_2_weight, proc_id_weight = 1, 1, 1, 3
ps_def, h_1, h_2 = 1, 1, 1


def _read_cms_angular_table(path):
    """
    Read a HEPData CMS 8TeV angular-coefficient file. It packs several
    A_i-vs-qT blocks (A0, A1, A2, A3, A4, A0-A2) separated by blank lines,
    each with its own header row "qT [GeV],qT [GeV] LOW,...,<label>,stat +,
    stat -,sys +,sys -". Returns {label: rows}, rows = list of float lists
    [qT_avg, qT_low, qT_high, value, stat+, stat-, sys+, sys-].
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    blocks = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("qT [GeV]"):
            label = line.split(",")[3]
            rows = []
            i += 1
            while i < len(lines) and lines[i]:
                rows.append([float(v) for v in lines[i].split(",")])
                i += 1
            blocks[label] = rows
        else:
            i += 1
    return blocks


def _add_cms_angular_points(ds, rows, y_min, y_max, proc_id):
    for i, row in enumerate(rows):
        qT_min, qT_max = row[1], row[2]
        xSec = row[3]

        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            ps_def_weight=ps_def_weight, h_1_weight=h_1_weight,
            h_2_weight=h_2_weight, proc_id_weight=proc_id_weight,
            s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            thFactor=1.,
            includeCuts=includeCuts,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            xSec=xSec,
            uncorrErr_0=(row[4] - row[5]) / 2.,
            uncorrErr_1=(row[6] - row[7]) / 2.,
        ))


blocks_00y10 = _read_cms_angular_table(path_to_data + "y0to1.csv")
blocks_10y21 = _read_cms_angular_table(path_to_data + "y1to21.csv")

#%% -- A0, 0<|y|<1
ds = DataSet.empty("DY", name="CMS8-A0-00y10", comment="CMS 8TeV angular coefficient A0, 0<|y|<1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_00y10["A0"], 0., 1., proc_id=20)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A0, 1<|y|<2.1
ds = DataSet.empty("DY", name="CMS8-A0-10y21", comment="CMS 8TeV angular coefficient A0, 1<|y|<2.1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_10y21["A0"], 1., 2.1, proc_id=20)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A1, 0<|y|<1
ds = DataSet.empty("DY", name="CMS8-A1-00y10", comment="CMS 8TeV angular coefficient A1, 0<|y|<1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_00y10["A1"], 0., 1., proc_id=21)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A1, 1<|y|<2.1
ds = DataSet.empty("DY", name="CMS8-A1-10y21", comment="CMS 8TeV angular coefficient A1, 1<|y|<2.1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_10y21["A1"], 1., 2.1, proc_id=21)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A2, 0<|y|<1
ds = DataSet.empty("DY", name="CMS8-A2-00y10", comment="CMS 8TeV angular coefficient A2, 0<|y|<1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_00y10["A2"], 0., 1., proc_id=22)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A2, 1<|y|<2.1
ds = DataSet.empty("DY", name="CMS8-A2-10y21", comment="CMS 8TeV angular coefficient A2, 1<|y|<2.1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_10y21["A2"], 1., 2.1, proc_id=22)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A3, 0<|y|<1
ds = DataSet.empty("DY", name="CMS8-A3-00y10", comment="CMS 8TeV angular coefficient A3, 0<|y|<1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_00y10["A3"], 0., 1., proc_id=23)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A3, 1<|y|<2.1
ds = DataSet.empty("DY", name="CMS8-A3-10y21", comment="CMS 8TeV angular coefficient A3, 1<|y|<2.1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_10y21["A3"], 1., 2.1, proc_id=23)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A4, 0<|y|<1
ds = DataSet.empty("DY", name="CMS8-A4-00y10", comment="CMS 8TeV angular coefficient A4, 0<|y|<1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_00y10["A4"], 0., 1., proc_id=24)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A4, 1<|y|<2.1
ds = DataSet.empty("DY", name="CMS8-A4-10y21", comment="CMS 8TeV angular coefficient A4, 1<|y|<2.1",
                    reference="1504.03512", normErr=[], isNormalized=False)
_add_cms_angular_points(ds, blocks_10y21["A4"], 1., 2.1, proc_id=24)
ds.save_csv(path_to_save + ds.name + ".csv")
