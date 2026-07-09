"""
Parsing of CMS 13 TeV Z/gamma* -> l+l- transverse-momentum spectrum, 5 rapidity
bins, absolute and normalized-to-1 versions.

Source: /data/arTeMiDe_Repository/data/CMS/CMS13-ydiff-abs-born.csv, CMS13-ydiff-norm-born.csv
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "CMS 13 y-diff" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(CMS13).py -- definitions
(process code, s, Q-range, y-ranges, thFactor, cuts, error split, the
hardcoded ptBins/FindBin bin-restoration table) are kept identical to the old
parsing.

NOTE: the raw tables only give each point's pT *bin-center label* (e.g. "1.5",
"15.0", ...), not the bin edges -- old parsing restored the (non-uniform)
edges from a hardcoded ptBins table (said to be "restored from the covariance
matrix file"). Ported verbatim below as _PT_BINS/_find_bin; not re-derived.

Datasets (arXiv:1909.04133):
    CMS13-00y04, CMS13-04y08, CMS13-08y12, CMS13-12y16, CMS13-16y24
        -- absolute cross section
    CMS13-00y04-norm, CMS13-04y08-norm, CMS13-08y12-norm, CMS13-12y16-norm, CMS13-16y24-norm
        -- normalized to 1
"""

import sys
import numpy as np
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.1876  # this script's own value (more precise than the 91. used elsewhere), ported as coded


### pT bins restored from the covariance matrix file (as in old parsing) -- the
### raw tables only give a bin-center label, not edges, and the bins are
### non-uniform, so this table cannot be re-derived from the data alone.
_PT_BINS = [
    [0., 1.], [1., 2.], [2., 3.], [3., 4.], [4., 5.], [5., 6.], [6., 7.], [7., 8.],
    [8., 9.], [9., 10.], [10., 11.], [11., 12.], [12., 13.], [13., 14.], [14., 16.],
    [16., 18.], [18., 20.], [20., 22.], [22., 25.], [25., 28.], [28., 32.], [32., 37.],
    [37., 43.], [43., 52.], [52., 65.], [65., 85.], [85., 120.], [120., 160.],
    [160., 190.], [190., 220.], [220., 250.], [250., 300.], [300., 400.], [400., 1500.],
]


def _find_bin(pt_center):
    for b in _PT_BINS:
        if b[0] <= pt_center < b[1]:
            return b
    raise ValueError(f"pT bin not found for center={pt_center}")


def _read_ydiff_blocks(path):
    """
    Read the 5 comma-separated blocks (one per rapidity bin) out of a
    CMS13-ydiff-*.csv file. Each block is a "$p_T$ [GeV], value, +, -" header
    line, 33 data rows, and a blank line; blocks start at (0-indexed) lines
    6, 41, 76, 111, 146. Returns a list of 5 lists of [pt_center, value, errp, errm] rows.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    starts = [5, 40, 75, 110, 145]  # 0-indexed header-line positions
    return [[[float(v) for v in lines[s + 1 + i].split(",")] for i in range(33)] for s in starts]


### given in the text (sec.6 / table 2)
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
s = 13000.**2
Q_min, Q_max = M_Z - 15., M_Z + 15.
includeCuts = True
cutParams = [25., 25., -2.4, 2.4]
### 2.5% luminosity uncertainty (table 2)
lumUncertainty = 0.025
### background estimations taken as fully correlated (table 2)
corrSys = 0.001

abs_blocks = _read_ydiff_blocks(path_to_data + "CMS/CMS13-ydiff-abs-born.csv")
norm_blocks = _read_ydiff_blocks(path_to_data + "CMS/CMS13-ydiff-norm-born.csv")

#%% -- 0.0 < |y| < 0.4, absolute
y_min, y_max = 0., 0.4
ds = DataSet.empty("DY", name="CMS13-00y04", comment="CMS 13TeV 0.0<|y|<0.4 absolute",
                    reference="arXiv:1909.04133", normErr=[lumUncertainty], isNormalized=False)

for i, (pt_center, xSec, errp, errm) in enumerate(abs_blocks[0]):
    qT_min, qT_max = _find_bin(pt_center)
    ### total quoted uncertainty includes the luminosity uncertainty; old parsing
    ### subtracts it in quadrature since normErr carries it separately
    unc = (errp - errm) / 2.
    uncorrErr_0 = np.sqrt(unc**2 - (xSec * lumUncertainty)**2)

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-00y04.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size; factor 2 to symmetrize y (as in old parsing)
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.4 < |y| < 0.8, absolute
y_min, y_max = 0.4, 0.8
ds = DataSet.empty("DY", name="CMS13-04y08", comment="CMS 13TeV 0.4<|y|<0.8 absolute",
                    reference="arXiv:1909.04133", normErr=[lumUncertainty], isNormalized=False)

for i, (pt_center, xSec, errp, errm) in enumerate(abs_blocks[1]):
    qT_min, qT_max = _find_bin(pt_center)
    unc = (errp - errm) / 2.
    uncorrErr_0 = np.sqrt(unc**2 - (xSec * lumUncertainty)**2)

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-04y08.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.8 < |y| < 1.2, absolute
y_min, y_max = 0.8, 1.2
ds = DataSet.empty("DY", name="CMS13-08y12", comment="CMS 13TeV 0.8<|y|<1.2 absolute",
                    reference="arXiv:1909.04133", normErr=[lumUncertainty], isNormalized=False)

for i, (pt_center, xSec, errp, errm) in enumerate(abs_blocks[2]):
    qT_min, qT_max = _find_bin(pt_center)
    unc = (errp - errm) / 2.
    uncorrErr_0 = np.sqrt(unc**2 - (xSec * lumUncertainty)**2)

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-08y12.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 1.2 < |y| < 1.6, absolute
y_min, y_max = 1.2, 1.6
ds = DataSet.empty("DY", name="CMS13-12y16", comment="CMS 13TeV 1.2<|y|<1.6 absolute",
                    reference="arXiv:1909.04133", normErr=[lumUncertainty], isNormalized=False)

for i, (pt_center, xSec, errp, errm) in enumerate(abs_blocks[3]):
    qT_min, qT_max = _find_bin(pt_center)
    unc = (errp - errm) / 2.
    uncorrErr_0 = np.sqrt(unc**2 - (xSec * lumUncertainty)**2)

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-12y16.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 1.6 < |y| < 2.4, absolute
y_min, y_max = 1.6, 2.4
ds = DataSet.empty("DY", name="CMS13-16y24", comment="CMS 13TeV 1.6<|y|<2.4 absolute",
                    reference="arXiv:1909.04133", normErr=[lumUncertainty], isNormalized=False)

for i, (pt_center, xSec, errp, errm) in enumerate(abs_blocks[4]):
    qT_min, qT_max = _find_bin(pt_center)
    unc = (errp - errm) / 2.
    uncorrErr_0 = np.sqrt(unc**2 - (xSec * lumUncertainty)**2)

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-16y24.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.0 < |y| < 0.4, normalized to 1
y_min, y_max = 0., 0.4
ds = DataSet.empty("DY", name="CMS13-00y04-norm", comment="CMS 13TeV 0.0<|y|<0.4 normalized",
                    reference="arXiv:1909.04133", normErr=[], isNormalized=True)

for i, (pt_center, xSec, errp, errm) in enumerate(norm_blocks[0]):
    qT_min, qT_max = _find_bin(pt_center)
    ### normalized data has no luminosity dependence to subtract (as in old parsing)
    uncorrErr_0 = (errp - errm) / 2.

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-00y04-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.4 < |y| < 0.8, normalized to 1
y_min, y_max = 0.4, 0.8
ds = DataSet.empty("DY", name="CMS13-04y08-norm", comment="CMS 13TeV 0.4<|y|<0.8 normalized",
                    reference="arXiv:1909.04133", normErr=[], isNormalized=True)

for i, (pt_center, xSec, errp, errm) in enumerate(norm_blocks[1]):
    qT_min, qT_max = _find_bin(pt_center)
    uncorrErr_0 = (errp - errm) / 2.

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-04y08-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 0.8 < |y| < 1.2, normalized to 1
y_min, y_max = 0.8, 1.2
ds = DataSet.empty("DY", name="CMS13-08y12-norm", comment="CMS 13TeV 0.8<|y|<1.2 normalized",
                    reference="arXiv:1909.04133", normErr=[], isNormalized=True)

for i, (pt_center, xSec, errp, errm) in enumerate(norm_blocks[2]):
    qT_min, qT_max = _find_bin(pt_center)
    uncorrErr_0 = (errp - errm) / 2.

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-08y12-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 1.2 < |y| < 1.6, normalized to 1
y_min, y_max = 1.2, 1.6
ds = DataSet.empty("DY", name="CMS13-12y16-norm", comment="CMS 13TeV 1.2<|y|<1.6 normalized",
                    reference="arXiv:1909.04133", normErr=[], isNormalized=True)

for i, (pt_center, xSec, errp, errm) in enumerate(norm_blocks[3]):
    qT_min, qT_max = _find_bin(pt_center)
    uncorrErr_0 = (errp - errm) / 2.

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-12y16-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 1.6 < |y| < 2.4, normalized to 1
y_min, y_max = 1.6, 2.4
ds = DataSet.empty("DY", name="CMS13-16y24-norm", comment="CMS 13TeV 1.6<|y|<2.4 normalized",
                    reference="arXiv:1909.04133", normErr=[], isNormalized=True)

for i, (pt_center, xSec, errp, errm) in enumerate(norm_blocks[4]):
    qT_min, qT_max = _find_bin(pt_center)
    uncorrErr_0 = (errp - errm) / 2.

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13-16y24-norm.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=2. * atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
        corrErr_0=xSec * corrSys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
