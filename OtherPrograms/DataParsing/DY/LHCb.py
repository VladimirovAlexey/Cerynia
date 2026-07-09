"""
Parsing of LHCb Z/gamma* -> mu+mu- transverse-momentum spectrum, 7/8/13 TeV.

Source: /data/arTeMiDe_Repository/data/LHCb/LHCb_7.dat, LHCb_8.dat, LHCb_13.dat,
        LHCb_13[2021].dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "LHCb" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py, plus the
y-integrated block of ReadDYdataFiles(LHCb13).py -- definitions (process
code, s, Q-range, y-range, thFactor, cuts, error split) are kept identical to
the old parsing.

Datasets:
    LHCb7        -- LHCb 7 TeV,  arXiv:1505.07024 (HEPData tab table)
    LHCb8        -- LHCb 8 TeV,  arXiv:1511.08039 (HEPData tab table)
    LHCb13(2016) -- LHCb 13 TeV, arXiv:1607.06495 (hand-typed table: no trailing
                    blank line, single symmetric error values instead of signed
                    dy+/dy- pairs). Renamed from the old parsing's plain "LHCb13"
                    to "LHCb13(2016)" per user instruction, to disambiguate from
                    the newer LHCb13(2021) below.
    LHCb13(2021) -- LHCb 13 TeV, 2021 update, arXiv:2112.07458, y-integrated
                    (2.0<y<4.5). Companion to the double-differential
                    LHCb13_dy.py (same paper); dataset name keeps the "(2021)"
                    label per user instruction (unlike LHCb13_dy, which had it
                    dropped).
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.  # forces <Q> = Z-boson mass rather than the (very close) bin midpoint,
           # matching the old parsing's explicit p["<Q>"]=M_Z override


def _read_hepdata_table(path, n_header, drop_trailing_blank=True):
    """
    Read a tab-separated HEPData .dat table.
    Each data row is "\\t x xlow xhigh y dy+ dy- [dy+ dy- ...]" -- the leading
    tab produces an empty first field, which is dropped. Returns a list of
    float rows [x, xlow, xhigh, y, dy+, dy-, ...].
    drop_trailing_blank=False for hand-typed files with no trailing blank line
    (e.g. LHCb_13.dat).
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    lines = lines[n_header:-1] if drop_trailing_blank else lines[n_header:]
    return [[float(v) for v in line.split("\t")[1:]] for line in lines]


#%%
# ============================================================================
# LHCb 7 TeV -- arXiv:1505.07024
# ============================================================================
rows = _read_hepdata_table(path_to_data + "LHCb/LHCb_7.dat", n_header=11)

### M(mu mu) window and PT(mu)>20 GeV cut, stated in the file header
Q_min, Q_max = 60., 120.
### 7 TeV
s = 7000.**2
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### ETARAP(MU) : 2.0 - 4.5, stated in the file header
y_min, y_max = 2., 4.5
includeCuts = True
cutParams = [20., 20., 2., 4.5]
### no luminosity/normalization error quoted (old parsing carried an explicit 0.0 entry)
normErr = []

### First (sigma_Z) block has 4 dy+/dy- pairs: stat, syst, beam, lumi (reduced
### offset 3). Old parsing kept stat+syst uncorrelated, beam+lumi correlated.
### Second block (f_FSR correction factor, offset 12) is unused, as in old parsing.

ds = DataSet.empty("DY", name="LHCb7", comment="LHCb 7TeV", reference="arXiv:1505.07024",
                    normErr=normErr, isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[3]
    dyp1, dym1, dyp2, dym2, dyp3, dym3, dyp4, dym4 = row[4:12]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"LHCb7.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### this data is not weighted by bin size (as in old parsing)
        thFactor=atmdeFactor,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        ### stat + syst, uncorrelated (as in old parsing)
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        ### beam + luminosity, correlated (as in old parsing)
        corrErr_0=(dyp3 - dym3) / 2.,
        corrErr_1=(dyp4 - dym4) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# LHCb 8 TeV -- arXiv:1511.08039
# ============================================================================
rows = _read_hepdata_table(path_to_data + "LHCb/LHCb_8.dat", n_header=11)

Q_min, Q_max = 60., 120.
### 8 TeV
s = 8000.**2
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
y_min, y_max = 2., 4.5
includeCuts = True
cutParams = [20., 20., 2., 4.5]
### no luminosity/normalization error quoted (old parsing carried an explicit 0.0 entry)
normErr = []

ds = DataSet.empty("DY", name="LHCb8", comment="LHCb 8TeV", reference="arXiv:1511.08039",
                    normErr=normErr, isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[3]
    dyp1, dym1, dyp2, dym2, dyp3, dym3, dyp4, dym4 = row[4:12]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"LHCb8.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### this data is not weighted by bin size (as in old parsing)
        thFactor=atmdeFactor,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        ### stat + syst, uncorrelated (as in old parsing)
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        ### beam + luminosity, correlated (as in old parsing)
        corrErr_0=(dyp3 - dym3) / 2.,
        corrErr_1=(dyp4 - dym4) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# LHCb 13 TeV, 2016 -- arXiv:1607.06495 (hand-typed, pending official HepData release)
# ============================================================================
rows = _read_hepdata_table(path_to_data + "LHCb/LHCb_13.dat", n_header=12, drop_trailing_blank=False)

Q_min, Q_max = 60., 120.
### 13 TeV
s = 13000.**2
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
y_min, y_max = 2., 4.5
includeCuts = True
cutParams = [20., 20., 2., 4.5]
### no luminosity/normalization error quoted (old parsing carried an explicit 0.0 entry)
normErr = []

### "(2016)" label added per user instruction, to disambiguate from LHCb13(2021) below
ds = DataSet.empty("DY", name="LHCb13(2016)", comment="LHCb 13TeV", reference="arXiv:1607.06495",
                    normErr=normErr, isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[1], row[2]
    xSec = row[3]
    ### this table gives single symmetric error values, not signed dy+/dy- pairs
    err1, err2, err3 = row[4], row[5], row[6]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"LHCb13(2016).{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=err1,
        uncorrErr_1=err2,
        corrErr_0=err3,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# LHCb 13 TeV, 2021 update, y-integrated -- arXiv:2112.07458
# ============================================================================
### hand-typed, comma-separated (NOT the tab-separated HEPData format the
### other LHCb files use) -- "ptMin , pTMax, xSec , stat. , sys. , lumi(2%)"
with open(path_to_data + "LHCb/LHCb_13[2021].dat") as f:
    _lines = [line.rstrip("\n") for line in f]
rows = [[float(v) for v in line.split(",")] for line in _lines[7:]]

Q_min, Q_max = 60., 120.
s = 13000.**2
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### ETARAP(MU) : 2.0 - 4.5, stated in the file header
y_min, y_max = 2., 4.5
includeCuts = True
cutParams = [20., 20., 2., 4.5]
### 2% luminosity uncertainty (table 1)
normErr = [0.02]

ds = DataSet.empty("DY", name="LHCb13(2021)", comment="LHCb 13TeV 2021 update",
                    reference="arXiv:2112.07458", normErr=normErr, isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[0], row[1]
    xSec = row[2]
    ### stat + syst, both treated uncorrelated (syst lightly correlated, ignored as in old parsing);
    ### lumi column (row[5]) is redundant with the global 2% normErr and unused per-point
    stat, sys = row[3], row[4]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"LHCb13(2021).{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=stat,
        uncorrErr_1=sys,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
