"""
Parsing of Tevatron Z/gamma* -> l+l- transverse-momentum spectra (CDF/D0, runs 1 & 2).

Source: /data/arTeMiDe_Repository/data/CDF_D0/*.dat  (HEPData tab-separated tables)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/*.csv

Ported from the "Tevatron" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q-range, y-range, thFactor, cuts, error split) are kept
identical to the old parsing.

Datasets:
    CDF1  -- CDF run 1, hep-ex/0001021
    CDF2  -- CDF run 2, arXiv:1207.7138
    D01   -- D0 run 1,  hep-ex/9907009
    D02   -- D0 run 2,  arXiv:0712.0803
    D02m  -- D0 run 2, muon channel, arXiv:1006.0618
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
# CDF run 1 -- hep-ex/0001021
# ============================================================================
rows = _read_hepdata_table(path_to_data + "CDF_D0/CDF_run1.dat", n_header=11)

###

### 3.9% luminocity error, as stated in the caption of table
normErr=[0.039]
### 1.8TeV
s= 1800.**2
### 2 page, left column: Only ee pairs in the mass range 66-116 GeV are used.
Q_min, Q_max  = 66., 116.
### The measurement is in full y- range. Maximum possible y=6.6... .I set +-10 for safety
y_min, y_max  = -10., 10.
### Process is Z-DY for p-pbar
ps_def,h_1, h_2,proc_id= 1,1,-1,3

ds = DataSet.empty("DY", name="CDF1", comment="CDF run1", reference="hep-ex/0001021",
                    normErr=normErr, isNormalized=False)


for i, (x, qT_min, qT_max, xSec, dyp, dym) in enumerate(rows):
    
    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor=(Q_max - Q_min)*(y_max-y_min)*(qT_max-qT_min)
    
    ds.add_point(dict(
        id=f"CDF1.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor/(qT_max-qT_min), ### the measurment in avaraged over pT 
        includeCuts=False,
        xSec=xSec, uncorrErr_0=(dyp - dym) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# CDF run 2 -- arXiv:1207.7138
# ============================================================================
rows = _read_hepdata_table(path_to_data + "CDF_D0/CDF_run2.dat", n_header=9)

### 5.8% luminosity error (see page 14, 2 line)
normErr = [0.058]
### 1.96 TeV
s = 1960.**2
### M(ee) window quoted in the HEPData record: 66-116 GeV
Q_min, Q_max = 66., 116.
### Full y-range measurement. Kinematic max |y|-range ~6.8 (2*ln(sqrt(s)/Q_min)). Set +-10 for safety
y_min, y_max = -10., 10.
### Process is Z-DY for p-pbar
ps_def, h_1, h_2, proc_id = 1, 1, -1, 3

ds = DataSet.empty("DY", name="CDF2", comment="CDF run2", reference="arXiv:1207.7138",
                    normErr=normErr, isNormalized=False)

for i, (x, qT_min, qT_max, xSec, dyp1, dym1, dyp2, dym2) in enumerate(rows):
    
    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor=(Q_max - Q_min)*(y_max-y_min)*(qT_max-qT_min)
    
    ds.add_point(dict(
        id=f"CDF2.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor/(qT_max-qT_min), ### the measurment in avaraged over pT 
        includeCuts=False,
        xSec=xSec,
        ### first published error band, uncorrelated (statistical)
        uncorrErr_0=(dyp1 - dym1) / 2.,
        ### second published error band, correlated systematic (as in old parsing)
        corrErr_0=(dyp2 - dym2) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# D0 run 1 -- hep-ex/9907009
# ============================================================================
rows = _read_hepdata_table(path_to_data + "CDF_D0/D0_run1.dat", n_header=9)

### 4.4% luminosity error (see Table 4 in paper)
normErr = [0.044]
### 1.8 TeV
s = 1800.**2
### Mass window given in paper
Q_min, Q_max = 75., 105.
### Full y-range measurement. Kinematic max |y|-range ~6.4. Set +-10 for safety
y_min, y_max = -10., 10.
### Process is Z-DY for p-pbar
ps_def, h_1, h_2, proc_id = 1, 1, -1, 3

ds = DataSet.empty("DY", name="D01", comment="D0 run1", reference="hep-ex/9907009",
                    normErr=normErr, isNormalized=False)

for i, (x, qT_min, qT_max, xSec, dyp, dym) in enumerate(rows):
    
    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor=(Q_max - Q_min)*(y_max-y_min)*(qT_max-qT_min)
    
    ds.add_point(dict(
        id=f"D01.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor/(qT_max-qT_min), ### the measurment in avaraged over pT 
        includeCuts=False,
        xSec=xSec, uncorrErr_0=(dyp - dym) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# D0 run 2 -- arXiv:0712.0803
# ============================================================================
rows = _read_hepdata_table(path_to_data + "CDF_D0/D0_run2.dat", n_header=9)

### Normalized (1/sig)dsig/dpT measurement -- no luminosity/normalization error
normErr = []
### 1.96 TeV
s = 1960.**2
### Mass window (paper-specified: invariant mass 70 < M (ee) < 110 GeV)
Q_min, Q_max = 70., 110.
### Full y-range measurement. Kinematic max |y|-range ~6.7 (2*ln(sqrt(s)/Q_min)). Set +-10 for safety
y_min, y_max = -10., 10.
### Process is Z-DY for p-pbar
ps_def, h_1, h_2, proc_id = 1, 1, -1, 3

ds = DataSet.empty("DY", name="D02", comment="D0 run2", reference="arXiv:0712.0803",
                    normErr=normErr, isNormalized=True)

for i, (x, qT_min, qT_max, xSec, dyp1, dym1, dyp2, dym2) in enumerate(rows):
    
    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor=(Q_max - Q_min)*(y_max-y_min)*(qT_max-qT_min)
    
    ds.add_point(dict(
        id=f"D02.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor/(qT_max-qT_min), ### the measurment in avaraged over pT 
        includeCuts=False,
        xSec=xSec,
        ### both published error bands are treated as uncorrelated (as in old parsing)
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%%
# ============================================================================
# D0 run 2, muon channel -- arXiv:1006.0618
# ============================================================================
rows = _read_hepdata_table(path_to_data + "CDF_D0/D0_run2m.dat", n_header=8)

### Normalized (1/sig)dsig/dpT measurement -- no luminosity/normalization error
normErr = []
### 1.96 TeV
s = 1960.**2
### Mass window (a dimuon mass in the range 65 < Mμμ < 115 GeV)
Q_min, Q_max = 65., 115.
### Full y-range measurement. Kinematic max |y|-range. Set +-10 for safety
y_min, y_max = -10., 10.
### Process is Z-DY for p-pbar
ps_def, h_1, h_2, proc_id = 1, 1, -1, 3
### Muon-channel fiducial cuts: p_T(mu) > 15 GeV, |eta(mu)| < 1.7 (see paper page 5)
includeCuts = True
cutParams = [15., 15., -1.7, 1.7]

ds = DataSet.empty("DY", name="D02m", comment="D0 run2 data for muons", reference="arXiv:1006.0618",
                    normErr=normErr, isNormalized=True)

for i, row in enumerate(rows):
    x, qT_min, qT_max, xSec = row[0], row[1], row[2], row[3]
    dyp1, dym1, dyp2, dym2, dyp3, dym3, dyp4, dym4, dyp5, dym5 = row[4:14]
    
    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor=(Q_max - Q_min)*(y_max-y_min)*(qT_max-qT_min)
    
    ds.add_point(dict(
        id=f"D02m.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor/(qT_max-qT_min), ### the measurment in avaraged over pT 
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        ### two uncorrelated error bands 
        uncorrErr_0=(dyp1 - dym1) / 2.,
        uncorrErr_1=(dyp2 - dym2) / 2.,
        ### three correlated (shape) systematic bands
        ### signs genuinely flip bin-to-bin, verified against the raw HEPData table
        corrErr_0=(dyp3 - dym3) / 2.,
        corrErr_1=(dyp4 - dym4) / 2.,
        corrErr_2=(dyp5 - dym5) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
