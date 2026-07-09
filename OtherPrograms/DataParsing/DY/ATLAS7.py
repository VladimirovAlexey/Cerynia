"""
Parsing of ATLAS 7 TeV Z/gamma* -> l+l- transverse-momentum spectrum, 3 rapidity bins.

Source: /data/arTeMiDe_Repository/data/ATLAS/ATLAS_7TeV_complete.dat  (HEPData tab-separated table)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "ATLAS 7" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q-range, y-ranges, thFactor, cuts, error split) are kept
identical to the old parsing.

Datasets (all from Table 3 of HepData record 8608, arXiv:1406.3660):
    A7-00y10  -- 0.0 < |y| < 1.0
    A7-10y20  -- 1.0 < |y| < 2.0
    A7-20y24  -- 2.0 < |y| < 2.4
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.  # forces <Q> = Z-boson mass rather than the (very close) bin midpoint,
           # matching the old parsing's explicit p["<Q>"]=M_Z override


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


#%%
# ============================================================================
# ATLAS 7 TeV, Table 3 -- arXiv:1406.3660
# ============================================================================
### Table 3 header runs to line 81 (three rapidity bins, Born+Dressed pair each)
rows = _read_hepdata_table(path_to_data + "ATLAS/ATLAS_7TeV_complete.dat", n_header=81)

### M(ll) window (stated in the begining of sec.7)
Q_min, Q_max = 66., 116.
### 7 TeV
s = 7000.**2
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### ATLAS 7 lepton fiducial cuts,  (stated in the begining of sec.7)
includeCuts = True
cutParams = [20., 20., -2.4, 2.4]
### Normalized (1/sig)dsig/dpT measurement -- no luminosity/normalization error
normErr = []

### Each rapidity bin occupies a (Born, Dressed) column pair (7 raw fields each);
### old parsing always used the first (Born) column of the pair -- raw offsets
### 3, 17, 31 into the reduced row (after dropping xdesc/x/xlow/xhigh).

#%% -- 0.0 < |y| < 1.0
y_min, y_max = 0., 1.
xSec_col = 3
ds = DataSet.empty("DY", name="A7-00y10", comment="ATLAS 7TeV 0.0<|y|<1.0", reference="arXiv:1406.3660",
                    normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A7-00y10.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
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

#%% -- 1.0 < |y| < 2.0
y_min, y_max = 1., 2.
xSec_col = 17
ds = DataSet.empty("DY", name="A7-10y20", comment="ATLAS 7TeV 1.0<|y|<2.0", reference="arXiv:1406.3660",
                    normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A7-10y20.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
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

#%% -- 2.0 < |y| < 2.4
y_min, y_max = 2., 2.4
xSec_col = 31
ds = DataSet.empty("DY", name="A7-20y24", comment="ATLAS 7TeV 2.0<|y|<2.4", reference="arXiv:1406.3660",
                    normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp1, dym1, dyp2, dym2, dyp3, dym3 = row[xSec_col + 1: xSec_col + 7]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"A7-20y24.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
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
