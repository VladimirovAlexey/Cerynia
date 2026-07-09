"""
Parsing of ATLAS 8 TeV Z/gamma* -> l+l- transverse-momentum spectrum,
normalized to 1/sigma (i.e. (1/sigma)*dsigma/dpT, not absolute cross section).

Same underlying measurement as ATLAS8.py (same 6 rapidity bins at the Z peak,
plus the 2 mass-sideband bins), but reading the "(1/sigma)dsigma/dpT" raw
tables instead of the "_abs" (absolute dsigma/dpT, in pb) ones.

Source: /data/arTeMiDe_Repository/data/ATLAS/ATLAS_8TeV_*.dat  (HEPData tables, HepData record 9030)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "ATLAS_norm" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(ATLAS_norm).py --
definitions (process code, s, Q-range, y-ranges, thFactor, cuts, error split)
are kept identical to the old parsing.

Datasets (arXiv:1512.02192):
    A8-00y04-norm   -- 66<Q<116 GeV, 0.0<|y|<0.4
    A8-04y08-norm   -- 66<Q<116 GeV, 0.4<|y|<0.8
    A8-08y12-norm   -- 66<Q<116 GeV, 0.8<|y|<1.2
    A8-12y16-norm   -- 66<Q<116 GeV, 1.2<|y|<1.6
    A8-16y20-norm   -- 66<Q<116 GeV, 1.6<|y|<2.0
    A8-20y24-norm   -- 66<Q<116 GeV, 2.0<|y|<2.4
    A8-46Q66-norm   -- 46<Q<66 GeV,  |y|<2.4 (rapidity-inclusive)
    A8-116Q150-norm -- 116<Q<150 GeV, |y|<2.4 (rapidity-inclusive)
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.  # forces <Q> = Z-boson mass rather than the (very close) bin midpoint,
           # matching the old parsing's explicit p["<Q>"]=M_Z override (Z-peak sets only)


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


### Same layout as the "_abs" files in this HepData record (9030): 6 columns
### per row (Electron dressed, Electron Born, Muon bare, Muon dressed, Muon
### Born, Combination Born), each a (y, dy+, dy-, dy+, dy-, dy+, dy-) block of
### 7 raw fields. Old parsing always used the last block ("Combination Born"),
### at reduced offset 38.
xSec_col = 38

### normalized (1/sigma)dsigma/dpT measurement -- no luminosity/normalization
### error (old parsing's normErr.append(lumUncertainty) line was commented out)
normErr = []
### 8 TeV
s = 8000.**2
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### ATLAS 8 lepton fiducial cuts (as in old parsing)
includeCuts = True
cutParams = [20., 20., -2.4, 2.4]

#%% -- 0.0 < |y| < 0.4, Z peak
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_66to116_00y04.dat", n_header=7)

Q_min, Q_max = 66., 116.
y_min, y_max = 0.0, 0.4

ds = DataSet.empty("DY", name="A8-00y04-norm", comment="ATLAS 8TeV 0.0<|y|<0.4 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-00y04-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.4 < |y| < 0.8, Z peak
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_66to116_04y08.dat", n_header=7)

Q_min, Q_max = 66., 116.
y_min, y_max = 0.4, 0.8

ds = DataSet.empty("DY", name="A8-04y08-norm", comment="ATLAS 8TeV 0.4<|y|<0.8 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-04y08-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.8 < |y| < 1.2, Z peak
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_66to116_08y12.dat", n_header=7)

Q_min, Q_max = 66., 116.
y_min, y_max = 0.8, 1.2

ds = DataSet.empty("DY", name="A8-08y12-norm", comment="ATLAS 8TeV 0.8<|y|<1.2 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-08y12-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 1.2 < |y| < 1.6, Z peak
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_66to116_12y16.dat", n_header=7)

Q_min, Q_max = 66., 116.
y_min, y_max = 1.2, 1.6

ds = DataSet.empty("DY", name="A8-12y16-norm", comment="ATLAS 8TeV 1.2<|y|<1.6 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-12y16-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 1.6 < |y| < 2.0, Z peak
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_66to116_16y20.dat", n_header=7)

Q_min, Q_max = 66., 116.
y_min, y_max = 1.6, 2.0

ds = DataSet.empty("DY", name="A8-16y20-norm", comment="ATLAS 8TeV 1.6<|y|<2.0 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-16y20-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 2.0 < |y| < 2.4, Z peak
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_66to116_20y24.dat", n_header=7)

Q_min, Q_max = 66., 116.
y_min, y_max = 2.0, 2.4

ds = DataSet.empty("DY", name="A8-20y24-norm", comment="ATLAS 8TeV 2.0<|y|<2.4 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-20y24-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 46 < Q < 66 GeV, rapidity-inclusive (|y|<2.4)
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_46to66.dat", n_header=7)

Q_min, Q_max = 46., 66.
y_min, y_max = -2.4, 2.4  ### already the full rapidity window; no symmetrize factor needed

ds = DataSet.empty("DY", name="A8-46Q66-norm", comment="ATLAS 8TeV 46<Q<66 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-46Q66-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,  ### off Z-peak: no Q_avg override, as in old parsing
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size, no symmetrize factor (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 116 < Q < 150 GeV, rapidity-inclusive (|y|<2.4)
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_8TeV_116to150.dat", n_header=7)

Q_min, Q_max = 116., 150.
y_min, y_max = -2.4, 2.4  ### already the full rapidity window; no symmetrize factor needed

ds = DataSet.empty("DY", name="A8-116Q150-norm", comment="ATLAS 8TeV 116<Q<150 normalized to 1/sigma",
                    reference="arXiv:1512.02192", normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A8-116Q150-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,  ### off Z-peak: no Q_avg override, as in old parsing
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size, no symmetrize factor (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        corrErr_0=(dyp3 - dym3) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
