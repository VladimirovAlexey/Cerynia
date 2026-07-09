"""
Parsing of COMPASS SIDIS multiplicities, deuteron target, isoscalar h+/h- final
states.

Source: /data/arTeMiDe_Repository/DataProcessor/DataLib/unpolSIDIS/Multiplicity_norm_byVMoos/compass(1706.01473).d.h{plus,minus}.datnew_fortran
        (already-preprocessed multiplicity tables, not raw COMPASS tables)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/SIDIS/*.csv

Ported from the "compass.deuteron.*" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSIDISdataFiles.py -- definitions
(process code, s, cuts, M_target/M_product, error split) are kept identical
to the old parsing. Process type is SIDIS -- see HERMES3D.py for the shared
schema notes (x/z/pT bins instead of qT/y).

NOTE: unlike HERMES, this raw table has no data-driven bin-average columns
(only bin edges), so Q_avg/x_avg/z_avg/pT_avg are left to Point.py's default
(bin midpoint) -- not set explicitly, matching old parsing (which never set
p["<Q>"] etc. for COMPASS).

NOTE: raw Q**2 and pT**2 columns are used directly in the old thFactor
(1/(pT_max**2-pT_min**2)/(z_max-z_min)/D11, D11 = data column 11, "xSec(DIS)")
-- structurally the same qT**2-Jacobian pattern seen in the DY branch's
E288/E605/E772. Per user confirmation, atmdeFactor = (Q_max-Q_min)*
(z_max-z_min)*(x_max-x_min)*(pT_max-pT_min) (using linear pT bounds) still
multiplies this, same standing rule as everywhere else, applied explicitly.

NOTE: old parsing had a special-cased fallback (thFactor=1) for xSec<1e-8,
same reasoning as HERMES. Audited both raw files (h+: 2332 rows, h-: 2332
rows): uniform 14 columns throughout, xSec is never <1e-8, D11 is never
zero -- the special case never actually triggers for COMPASS. No rows
dropped here (unlike HERMES); ported the single thFactor expression only.

NOTE: old parsing's `p["m_product"]=m_pion` uses a lowercase key ("m_product"
not "M_product") -- a likely typo that silently fell back to Point's default,
which is 0.139 (pion mass) anyway, so no numeric difference. Set M_product
explicitly and correctly-spelled here.

Datasets (1709.07374):
    compass.d.h+
    compass.d.h-
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/DataProcessor/DataLib/unpolSIDIS/Multiplicity_norm_byVMoos/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/SIDIS/"

M_proton = 0.938
m_pion = 0.139

### COMPASS beam energy 160 GeV on a fixed target -> s = 2*E_beam*M_target + M_target**2
s = 2 * 160. * M_proton + M_proton**2
### COMPASS fiducial cuts: [yMin, yMax, W2min, W2max] (as in old parsing)
cutParams = [0.1, 0.9, 25., 10000.]
### Process is SIDIS is 2001
### The deuteron corresponds to h=12
### the produced hadron is h=pi+K = 12 (for h+) or -12 (for h-)
ps_def, h_1, proc_id = 1, 12, 2001


def _read_compass_table(path, n_header=3):
    """
    Read a COMPASS SIDIS multiplicity table. Columns: x_min, x_max, Q2_min,
    Q2_max, z_min, z_max, pT2_min, pT2_max, mult, stat, relSyst, xSec(DIS),
    xSec(DIS-Err)+, xSec(DIS-Err)-. No leading row-index field (unlike HERMES).
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    for line in lines[n_header:]:
        parts = line.split()
        if not parts:
            continue
        rows.append([float(v) for v in parts])
    return rows


#%% -- deuteron target, h+
rows = _read_compass_table(path_to_data + "compass(1706.01473).d.hplus.datnew_fortran")

ds = DataSet.empty("SIDIS", name="compass.d.h+", comment="COMPASS isoscalar h+. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1709.07374", normErr=[], isNormalized=False)

for i, row in enumerate(rows):
    x_min, x_max = row[0], row[1]
    Q_min, Q_max = sqrt(row[2]), sqrt(row[3])
    z_min, z_max = row[4], row[5]
    pT_min, pT_max = sqrt(row[6]), sqrt(row[7])
    xSec = row[8]
    D11 = row[11]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (z_max - z_min) * (x_max - x_min) * (pT_max - pT_min)

    ds.add_point(dict(
        id=f"compass.d.h+.{i}", ps_def=ps_def, h_1=h_1, h_2=12, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, x_min=x_min, x_max=x_max, z_min=z_min, z_max=z_max,
        pT_min=pT_min, pT_max=pT_max,
        M_target=M_proton, M_product=m_pion,
        ### divide by bin size (pT in quadrature, matching the raw pT**2 columns), multiply
        ### by DIS xSec normalization (as in old parsing)
        thFactor=atmdeFactor / (pT_max**2 - pT_min**2) / (z_max - z_min) / D11,
        includeCuts=True,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=row[9],
        ### relative systematic in the raw file; converted to absolute (as in old parsing)
        uncorrErr_1=row[10] * xSec,
        ### DIS normalization uncertainty, treated as correlated (as in old parsing)
        corrErr_0=(row[12] + row[13]) / 2. / D11 * xSec,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, h-
rows = _read_compass_table(path_to_data + "compass(1706.01473).d.hminus.datnew_fortran")

ds = DataSet.empty("SIDIS", name="compass.d.h-", comment="COMPASS isoscalar h-. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1709.07374", normErr=[], isNormalized=False)

for i, row in enumerate(rows):
    x_min, x_max = row[0], row[1]
    Q_min, Q_max = sqrt(row[2]), sqrt(row[3])
    z_min, z_max = row[4], row[5]
    pT_min, pT_max = sqrt(row[6]), sqrt(row[7])
    xSec = row[8]
    D11 = row[11]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (z_max - z_min) * (x_max - x_min) * (pT_max - pT_min)

    ds.add_point(dict(
        id=f"compass.d.h-.{i}", ps_def=ps_def, h_1=h_1, h_2=-12, proc_id=proc_id,
        s=s, Q_min=Q_min, Q_max=Q_max, x_min=x_min, x_max=x_max, z_min=z_min, z_max=z_max,
        pT_min=pT_min, pT_max=pT_max,
        M_target=M_proton, M_product=m_pion,
        thFactor=atmdeFactor / (pT_max**2 - pT_min**2) / (z_max - z_min) / D11,
        includeCuts=True,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=row[9],
        uncorrErr_1=row[10] * xSec,
        corrErr_0=(row[12] + row[13]) / 2. / D11 * xSec,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
