"""
Parsing of CMS 13 TeV Z/gamma* -> l+l- transverse-momentum spectrum, 5 wide
Q-windows (preliminary result).

Source: /data/arTeMiDe_Repository/data/CMS/CMS13_50to76.dat, CMS13_76to106.dat,
        CMS13_106to170.dat, CMS13_170to350.dat, CMS13_350to1000.dat
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/CMS13_dQ_*.csv

Ported from the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(CMS13-largeQ).py --
definitions (process code, s, y-range, thFactor, cuts, error split, the
hardcoded per-point FSR correction factors) are kept identical to the old
parsing. Separate script from CMS13.py because the source data is completely
different (hand-typed preliminary CMS-PAS table vs. the published HepData
CMS13-ydiff-*.csv tables).

Raw format: hand-typed from CMS-PAS-SMP-20-003 (not a HEPData export),
comma-separated "ptMin, ptMax, xSec, uncer.total(%), <breakdown by source...>"
-- only the total percentage uncertainty (4th column) is used; the per-source
breakdown columns are unused, exactly as old parsing.

NOTE: the "fsr" correction factors below are hardcoded values from the old
script ("FSR factors computed by Louis Moureaux ... I extract them by taking
the ratio of nofsr/fsr") -- ported verbatim, not re-derived; there is no way
to recompute them from the raw data files alone.

Datasets (CMS-PAS-SMP-20-003):
    CMS13_dQ_50to76, CMS13_dQ_76to106, CMS13_dQ_106to170,
    CMS13_dQ_170to350, CMS13_dQ_350to1000
"""

import sys
import numpy as np
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"


def _read_cms13_largeq(path, n_header=10):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    return [[float(v) for v in line.split(",")] for line in lines[n_header:]]


### given in the text
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
s = 13000.**2
includeCuts = True
cutParams = [25., 20., -2.4, 2.4]
### 1.2% luminosity uncertainty
lumUncertainty = 0.012
y_min, y_max = -2.4, 2.4

#### FSR factors computed by Louis Moureaux (ratio of nofsr/fsr); ported verbatim
_FSR = {
    "50to76": [
        1.00098537, 0.97497963, 0.92974446, 0.87653573, 0.80923609, 0.73646756,
        0.66225815, 0.5801822, 0.51508677, 0.51669914, 0.59243914, 0.68497092,
        0.72740284, 0.84009018, 0.82205728, 0.9087,
    ],
    "76to106": [
        1.05285829, 1.05091075, 1.04498482, 1.03978998, 1.03490338, 1.02954884,
        1.02554925, 1.02210408, 1.02117097, 1.01910361, 1.01701754, 1.01836032,
        1.01750745, 1.01730813, 1.01867058, 1.02046079, 1.02091208, 1.02347862,
        1.02336165, 1.02557772, 1.0269174, 1.02751717, 1.0293184, 1.03087537,
        1.03446688, 1.0432738, 0.9932898, 1.025378, 1.02528262, 1.02589994,
        1.02191976, 1.02025379, 1.02884329, 1.02806738, 1.0233248, 1.01200851,
        1.02687311,
    ],
    "106to170": [
        1.07275747, 1.066437, 1.06056855, 1.05376852, 1.05392052, 1.04900427,
        1.04836114, 1.04849858, 1.04644792, 1.04296663, 1.04972161, 1.0434394,
        1.04831678, 1.68206705, 1.10846788, 1.04015155,
    ],
    "170to350": [
        1.08153606, 1.06157922, 1.06100099, 1.06240545, 1.039029, 1.04436067,
        1.03852346, 1.02823426, 1.02717946, 1.0260605, 1.0184391, 1.00729932,
        1.06753421,
    ],
    "350to1000": [
        1.07715219, 1.08831779, 1.03406113, 1.05885152, 1.09046837, 1.02711162,
        1.04009832, 1.03338694, 1.04638257, 1.01053939, 1.01408975, 1.00578608,
        1.04662236,
    ],
}

#%% -- 50 < Q < 76 GeV
Q_min, Q_max = 50., 76.
tag = "50to76"
rows = _read_cms13_largeq(path_to_data + f"CMS/CMS13_{tag}.dat")
ds = DataSet.empty("DY", name=f"CMS13_dQ_{tag}", comment="CMS 13TeV 2021 preliminary",
                    reference="CMS-PAS-SMP-20-003", normErr=[lumUncertainty], isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[0], row[1]
    xSec = row[2]
    pct_total = row[3]
    ### total % uncertainty includes the 1.2% luminosity uncertainty; old parsing
    ### subtracts it in quadrature since normErr carries it separately. Remaining
    ### breakdown-by-source columns are lightly correlated and ignored, as in old parsing.
    uncorrErr_0 = xSec * np.sqrt(pct_total**2 - 1.2**2) * 0.01

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13_dQ_{tag}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### divide by bin size and by the point's FSR correction factor (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min) / _FSR[tag][i],
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 76 < Q < 106 GeV
Q_min, Q_max = 76., 106.
tag = "76to106"
rows = _read_cms13_largeq(path_to_data + f"CMS/CMS13_{tag}.dat")
ds = DataSet.empty("DY", name=f"CMS13_dQ_{tag}", comment="CMS 13TeV 2021 preliminary",
                    reference="CMS-PAS-SMP-20-003", normErr=[lumUncertainty], isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[0], row[1]
    xSec = row[2]
    pct_total = row[3]
    uncorrErr_0 = xSec * np.sqrt(pct_total**2 - 1.2**2) * 0.01

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13_dQ_{tag}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor / (qT_max - qT_min) / _FSR[tag][i],
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 106 < Q < 170 GeV
Q_min, Q_max = 106., 170.
tag = "106to170"
rows = _read_cms13_largeq(path_to_data + f"CMS/CMS13_{tag}.dat")
ds = DataSet.empty("DY", name=f"CMS13_dQ_{tag}", comment="CMS 13TeV 2021 preliminary",
                    reference="CMS-PAS-SMP-20-003", normErr=[lumUncertainty], isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[0], row[1]
    xSec = row[2]
    pct_total = row[3]
    uncorrErr_0 = xSec * np.sqrt(pct_total**2 - 1.2**2) * 0.01

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13_dQ_{tag}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor / (qT_max - qT_min) / _FSR[tag][i],
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 170 < Q < 350 GeV
Q_min, Q_max = 170., 350.
tag = "170to350"
rows = _read_cms13_largeq(path_to_data + f"CMS/CMS13_{tag}.dat")
ds = DataSet.empty("DY", name=f"CMS13_dQ_{tag}", comment="CMS 13TeV 2021 preliminary",
                    reference="CMS-PAS-SMP-20-003", normErr=[lumUncertainty], isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[0], row[1]
    xSec = row[2]
    pct_total = row[3]
    uncorrErr_0 = xSec * np.sqrt(pct_total**2 - 1.2**2) * 0.01

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13_dQ_{tag}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor / (qT_max - qT_min) / _FSR[tag][i],
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- 350 < Q < 1000 GeV
Q_min, Q_max = 350., 1000.
tag = "350to1000"
rows = _read_cms13_largeq(path_to_data + f"CMS/CMS13_{tag}.dat")
ds = DataSet.empty("DY", name=f"CMS13_dQ_{tag}", comment="CMS 13TeV 2021 preliminary",
                    reference="CMS-PAS-SMP-20-003", normErr=[lumUncertainty], isNormalized=False)

for i, row in enumerate(rows):
    qT_min, qT_max = row[0], row[1]
    xSec = row[2]
    pct_total = row[3]
    uncorrErr_0 = xSec * np.sqrt(pct_total**2 - 1.2**2) * 0.01

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"CMS13_dQ_{tag}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        thFactor=atmdeFactor / (qT_max - qT_min) / _FSR[tag][i],
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=uncorrErr_0,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
