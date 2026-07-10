"""
Parsing of ATLAS 8 TeV Z-boson angular coefficients A0-A7 (Collins-Soper frame),
3 rapidity bins, near the Z peak.

Source: /data/arTeMiDe_Repository/data/ATLAS/Zangular_1606_00689/Table*.csv
        (HEPData record for arXiv:1606.00689, "regularised" tables)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY_angular/*.csv

Ported from the "A8_A*" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(angular).py --
definitions (process code, s, Q-range, y-ranges, thFactor, cuts, error split)
are kept identical to the old parsing.

NOTE: proc_id numbering (20-27 for A0-A7) confirmed correct by the user --
matches modern artemide's TMDX_DY process-code convention.

NOTE: the old library's "A8_Auu_*" datasets are NOT ported. They read the
exact same raw tables (and the exact same numeric column) as the A4 datasets
(Table36/56/75) -- i.e. they were byte-for-byte duplicates of A4 mislabeled
as the Auu/normalization observable, not real cross-section data. In the old
scheme these were used to build a weighting cross-section for the angular
coefficients; per user instruction, this is no longer ported as a separate
dataset. Instead, every point here carries the new weight-process columns
(ps_def_weight, h_1_weight, h_2_weight, proc_id_weight) = (1,1,1,3), telling
artemide which process to use as the normalization/weighting cross section
when computing the theoretical prediction for each ratio observable.

NOTE: thFactor=1 (ported verbatim, not given the usual atmdeFactor bin-width
correction) -- these are dimensionless ratio observables, and the weight-
process columns above let artemide bin-average numerator and denominator
consistently, so no extra Jacobian is applied here.

Datasets (arXiv:1606.00689), 22 total (A1 and A6 are missing at 2<|y|<3.5,
not measured in the source):
    A8-A0-00y10, A8-A0-10y20, A8-A0-20y35
    A8-A1-00y10, A8-A1-10y20
    A8-A2-00y10, A8-A2-10y20, A8-A2-20y35
    A8-A3-00y10, A8-A3-10y20, A8-A3-20y35
    A8-A4-00y10, A8-A4-10y20, A8-A4-20y35
    A8-A5-00y10, A8-A5-10y20, A8-A5-20y35
    A8-A6-00y10, A8-A6-10y20
    A8-A7-00y10, A8-A7-10y20, A8-A7-20y35
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/ATLAS/Zangular_1606_00689/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY_angular/"

M_Z = 91.1876  # forces <Q> = Z-boson mass, matching old parsing's explicit p["<Q>"]=M_Z

### 8 TeV, near the Z peak (as in old parsing)
s = 8000.**2
Q_min, Q_max = 80., 100.
### ATLAS 8 lepton fiducial cuts, unused since includeCuts=False (as in old parsing)
includeCuts = False
cutParams = [20., 20., 2., 4.5]
### weighting process: plain unpolarized DY cross section (ps_def,h_1,h_2,proc_id),
### used by artemide to build the normalization for these ratio observables
ps_def_weight, h_1_weight, h_2_weight, proc_id_weight = 1, 1, 1, 3
ps_def, h_1, h_2 = 1, 1, 1


def _read_atlas_angular_table(path):
    """
    Read a HEPData "regularised" A_i-vs-pT table. Columns: PT_avg, PT_low,
    PT_high, A_i, stat+, stat-, sys+, sys-, reg_bias+, reg_bias-.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    for line in lines:
        if not line or line[0] in ("#", "P"):
            continue
        rows.append([float(v) for v in line.split(",")])
    return rows


def _add_atlas_angular_points(ds, rows, y_min, y_max, proc_id):
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
            uncorrErr_2=(row[8] - row[9]) / 2.,
        ))


#%% -- A0, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table32.csv")
ds = DataSet.empty("DY", name="A8-A0-00y10", comment="ATLAS 8TeV angular coefficient A0, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=20)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A0, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table52.csv")
ds = DataSet.empty("DY", name="A8-A0-10y20", comment="ATLAS 8TeV angular coefficient A0, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=20)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A0, 2<|y|<3.5
rows = _read_atlas_angular_table(path_to_data + "Table74.csv")
ds = DataSet.empty("DY", name="A8-A0-20y35", comment="ATLAS 8TeV angular coefficient A0, 2<|y|<3.5",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 2., 3.5, proc_id=20)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A1, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table33.csv")
ds = DataSet.empty("DY", name="A8-A1-00y10", comment="ATLAS 8TeV angular coefficient A1, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=21)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A1, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table53.csv")
ds = DataSet.empty("DY", name="A8-A1-10y20", comment="ATLAS 8TeV angular coefficient A1, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=21)
ds.save_csv(path_to_save + ds.name + ".csv")

#################################################################
##### there is no A1 measurement for 2<|y|<3.5

#%% -- A2, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table34.csv")
ds = DataSet.empty("DY", name="A8-A2-00y10", comment="ATLAS 8TeV angular coefficient A2, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=22)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A2, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table54.csv")
ds = DataSet.empty("DY", name="A8-A2-10y20", comment="ATLAS 8TeV angular coefficient A2, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=22)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A2, 2<|y|<3.5
rows = _read_atlas_angular_table(path_to_data + "Table66.csv")
ds = DataSet.empty("DY", name="A8-A2-20y35", comment="ATLAS 8TeV angular coefficient A2, 2<|y|<3.5",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 2., 3.5, proc_id=22)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A3, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table35.csv")
ds = DataSet.empty("DY", name="A8-A3-00y10", comment="ATLAS 8TeV angular coefficient A3, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=23)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A3, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table55.csv")
ds = DataSet.empty("DY", name="A8-A3-10y20", comment="ATLAS 8TeV angular coefficient A3, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=23)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A3, 2<|y|<3.5
rows = _read_atlas_angular_table(path_to_data + "Table67.csv")
ds = DataSet.empty("DY", name="A8-A3-20y35", comment="ATLAS 8TeV angular coefficient A3, 2<|y|<3.5",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 2., 3.5, proc_id=23)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A4, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table36.csv")
ds = DataSet.empty("DY", name="A8-A4-00y10", comment="ATLAS 8TeV angular coefficient A4, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=24)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A4, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table56.csv")
ds = DataSet.empty("DY", name="A8-A4-10y20", comment="ATLAS 8TeV angular coefficient A4, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=24)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A4, 2<|y|<3.5
rows = _read_atlas_angular_table(path_to_data + "Table75.csv")
ds = DataSet.empty("DY", name="A8-A4-20y35", comment="ATLAS 8TeV angular coefficient A4, 2<|y|<3.5",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 2., 3.5, proc_id=24)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A5, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table37.csv")
ds = DataSet.empty("DY", name="A8-A5-00y10", comment="ATLAS 8TeV angular coefficient A5, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=25)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A5, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table57.csv")
ds = DataSet.empty("DY", name="A8-A5-10y20", comment="ATLAS 8TeV angular coefficient A5, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=25)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A5, 2<|y|<3.5
rows = _read_atlas_angular_table(path_to_data + "Table68.csv")
ds = DataSet.empty("DY", name="A8-A5-20y35", comment="ATLAS 8TeV angular coefficient A5, 2<|y|<3.5",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 2., 3.5, proc_id=25)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A6, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table38.csv")
ds = DataSet.empty("DY", name="A8-A6-00y10", comment="ATLAS 8TeV angular coefficient A6, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=26)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A6, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table58.csv")
ds = DataSet.empty("DY", name="A8-A6-10y20", comment="ATLAS 8TeV angular coefficient A6, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=26)
ds.save_csv(path_to_save + ds.name + ".csv")

#################################################################
##### there is no A6 measurement for 2<|y|<3.5

#%% -- A7, 0<|y|<1
rows = _read_atlas_angular_table(path_to_data + "Table39.csv")
ds = DataSet.empty("DY", name="A8-A7-00y10", comment="ATLAS 8TeV angular coefficient A7, 0<|y|<1",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 0., 1., proc_id=27)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A7, 1<|y|<2
rows = _read_atlas_angular_table(path_to_data + "Table59.csv")
ds = DataSet.empty("DY", name="A8-A7-10y20", comment="ATLAS 8TeV angular coefficient A7, 1<|y|<2",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 1., 2., proc_id=27)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- A7, 2<|y|<3.5
rows = _read_atlas_angular_table(path_to_data + "Table69.csv")
ds = DataSet.empty("DY", name="A8-A7-20y35", comment="ATLAS 8TeV angular coefficient A7, 2<|y|<3.5",
                    reference="1606.00689", normErr=[], isNormalized=False)
_add_atlas_angular_points(ds, rows, 2., 3.5, proc_id=27)
ds.save_csv(path_to_save + ds.name + ".csv")
