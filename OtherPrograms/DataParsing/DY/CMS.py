"""
Parsing of CMS Z/gamma* -> l+l- transverse-momentum spectrum, 7 and 8 TeV.

Source: /data/arTeMiDe_Repository/data/CMS/CMS_7.dat, CMS_8.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "CMS" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q-range, y-range, thFactor, cuts, error split) are kept
identical to the old parsing.

Datasets:
    CMS7 -- CMS 7 TeV, arXiv:1110.4973 (HEPData tab table)
    CMS8 -- CMS 8 TeV, arXiv:1606.05864 (YODA Scatter2D table -- different raw format)
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
# CMS 7 TeV -- arXiv:1110.4973
# ============================================================================
rows = _read_hepdata_table(path_to_data + "CMS/CMS_7.dat", n_header=11)

### M(ll) window, lepton |eta|<2.1, pT>20 GeV cuts -- all stated in the file header
Q_min, Q_max = 60., 120.
### 7 TeV
s = 7000.**2
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
y_min, y_max = -2.1, 2.1
includeCuts = True
cutParams = [20., 20., -2.1, 2.1]
### no luminosity/normalization error quoted (old parsing carried an explicit 0.0 entry)
normErr = []

### Header lists 3 channels (muon, electron, combined lepton) as 3 (y,dy+,dy-)
### blocks; old parsing used the third ("combined") block, at reduced offset 9.
xSec_col = 9

ds = DataSet.empty("DY", name="CMS7", comment="CMS 7TeV", reference="arXiv:1110.4973",
                    normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[xSec_col]
    dyp, dym = row[xSec_col + 1], row[xSec_col + 2]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS7.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(dyp - dym) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# CMS 8 TeV -- arXiv:1606.05864  (YODA Scatter2D table, different raw format)
# ============================================================================
### Columns are "xval xerr- xerr+ yval yerr- yerr+", no leading tab (unlike the
### HEPData xdesc format). Header is 6 lines (BEGIN/IsRef/Path/Title/Type/# xval);
### footer is "END YODA_SCATTER2D" + a trailing blank line (2 lines dropped).
### yerr-/yerr+ are both positive magnitudes here (not signed +/-).
with open(path_to_data + "CMS/CMS_8.dat") as f:
    _lines = [l.rstrip("\n") for l in f.readlines()]
rows = [[float(v) for v in l.split("\t")] for l in _lines[6:-2]]

### same analysis window/cuts as CMS7
Q_min, Q_max = 60., 120.
### 8 TeV
s = 8000.**2
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
y_min, y_max = -2.1, 2.1
includeCuts = True
cutParams = [20., 20., -2.1, 2.1]
### no luminosity/normalization error quoted (old parsing carried an explicit 0.0 entry)
normErr = []

ds = DataSet.empty("DY", name="CMS8", comment="CMS 8TeV", reference="arXiv:1606.05864",
                    normErr=normErr, isNormalized=True)

for i, (xval, xerrm, xerrp, yval, yerrm, yerrp) in enumerate(rows):
    qT_min, qT_max = xval - xerrm, xval + xerrp

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS8.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=yval,
        ### yerr-/yerr+ are both positive magnitudes in this format; average them
        uncorrErr_0=(yerrm + yerrp) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
