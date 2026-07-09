"""
Parsing of HERMES SIDIS multiplicities (zxpt-3D binning: z, x, pT), proton and
deuteron targets, pi+/pi-/K+/K- final states, with and without vector-meson
background subtraction (vmsub).

Source: /data/arTeMiDe_Repository/DataProcessor/DataLib/unpolSIDIS/Multiplicity_norm_byVMoos/hermes.*.datnew_fortran
        (already-preprocessed multiplicity tables, not raw HERMES tables)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/SIDIS/*.csv

Ported from the "hermes.*" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSIDISdataFiles.py -- definitions
(process code, s, cuts, M_target/M_product, error split) are kept identical
to the old parsing. Process type is SIDIS (not DY): points carry x/z/pT bins
instead of qT/y, per Cerynia's Point.py schema.

Datasets renamed per user instruction to match the COMPASS branch's naming
convention (see COMPASS.py): old "hermes.p.vmsub.zxpt.pi+" -> "hermes3D.p.pi+"
(vmsub is now the unlabeled/default case), old "hermes.p.no-vmsub.zxpt.pi+"
-> "hermes3D.p.pi+.no-vmsub" (no-vmsub gets an explicit suffix instead).

NOTE: thFactor is read from a column in the raw table (the "DIS normalization"
value, MSHT20nnlo-based), not purely computed from bin widths -- old formula
(HERMES): 1/(pT_max-pT_min)/(z_max-z_min)/D15, where D15 = data column 15
("xSec(DIS)"). Per user confirmation, atmdeFactor = (Q_max-Q_min)*(z_max-z_min)
*(x_max-x_min)*(pT_max-pT_min) still multiplies this (same standing rule as
DY), applied explicitly rather than pre-simplified.

NOTE: old parsing had a special-cased fallback (thFactor=1) for rows with
xSec<1e-8, to avoid a 0/0 (all raw fields, including D15, are exactly zero
for these rows). Audited every raw file: these rows are not real data --
they are entirely-zero placeholder rows for unpopulated (x,z,pT) grid cells,
also missing their very last column (17 fields instead of 18). Per user
confirmation, such rows are dropped entirely (not added as points), and the
remaining (real) rows all use the same thFactor expression -- no more
special-casing needed. Verified: every kept row has D15 != 0 and xSec >= 1e-8
in every one of the 16 raw files.

NOTE: M_target is M_proton (0.938) for BOTH proton and deuteron targets in
old parsing (deuteron target's per-nucleon SIDIS structure functions use the
single-nucleon mass) -- ported verbatim, not a bug.

Datasets (1212.5407), 16 total = 2 targets x 2 (vmsub/no-vmsub) x 4 hadrons:
    hermes3D.p.pi+, hermes3D.p.pi-, hermes3D.p.k+, hermes3D.p.k-
    hermes3D.p.pi+.no-vmsub, hermes3D.p.pi-.no-vmsub, hermes3D.p.k+.no-vmsub, hermes3D.p.k-.no-vmsub
    hermes3D.d.pi+, hermes3D.d.pi-, hermes3D.d.k+, hermes3D.d.k-
    hermes3D.d.pi+.no-vmsub, hermes3D.d.pi-.no-vmsub, hermes3D.d.k+.no-vmsub, hermes3D.d.k-.no-vmsub
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/DataProcessor/DataLib/unpolSIDIS/Multiplicity_norm_byVMoos/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/SIDIS/"

M_proton = 0.938
m_pion = 0.139
m_kaon = 0.494

### HERMES beam energy 27.6 GeV on a fixed target -> s = 2*E_beam*M_target + M_target**2
s = 2 * 27.6 * M_proton + M_proton**2
### HERMES fiducial cuts: [yMin, yMax, W2min, W2max] (as in old parsing)
cutParams = [0.1, 0.85, 10., 10000.]


def _read_hermes_table(path, n_header=5):
    """
    Read a HERMES zxpt-3D multiplicity table. Columns (after dropping the
    leading row-index field): mult, stat, syst, Q2_avg, Q2_min, Q2_max,
    x_avg, x_min, x_max, z_avg, z_min, z_max, pT_avg, pT_min, pT_max,
    xSec(DIS), xSec(DIS-Err)+, xSec(DIS-Err)-.
    Rows that are entirely-zero placeholders for unpopulated (x,z,pT) grid
    cells (missing their last column too -- 17 fields instead of 18) are
    dropped; see module docstring.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    for line in lines[n_header:]:
        parts = line.split()
        if not parts:
            continue
        vals = [float(v) for v in parts[1:]]
        if len(vals) < 18:
            continue  # all-zero placeholder row
        rows.append(vals)
    return rows


def _add_hermes_points(ds, rows, ps_def, h_1, h_2, proc_id, M_product):
    for i, row in enumerate(rows):
        xSec = row[0]
        Q_avg, Q_min, Q_max = sqrt(row[3]), sqrt(row[4]), sqrt(row[5])
        x_avg, x_min, x_max = row[6], row[7], row[8]
        z_avg, z_min, z_max = row[9], row[10], row[11]
        pT_avg, pT_min, pT_max = row[12], row[13], row[14]
        D15 = row[15]

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (z_max - z_min) * (x_max - x_min) * (pT_max - pT_min)

        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg,
            x_min=x_min, x_max=x_max, x_avg=x_avg,
            z_min=z_min, z_max=z_max, z_avg=z_avg,
            pT_min=pT_min, pT_max=pT_max, pT_avg=pT_avg,
            M_target=M_proton, M_product=M_product,
            ### divide by bin size, multiply by DIS xSec normalization (as in old parsing)
            thFactor=atmdeFactor / (pT_max - pT_min) / (z_max - z_min) / D15,
            includeCuts=True,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            xSec=xSec,
            uncorrErr_0=row[1],
            uncorrErr_1=row[2],
            ### DIS normalization uncertainty, treated as correlated (as in old parsing)
            corrErr_0=(row[16] + row[17]) / 2. / D15 * xSec,
        ))


### Process is SIDIS is 2001
### The proton target corresponds to h=1, deuteron to h=12
### the produced hadron is (pi+=1, pi-=-1, K+=2, K-=-2);
ps_def, proc_id = 1, 2001

#%% -- proton target, pi+, vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.vmsub.mults_piplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.pi+", comment="HERMES proton-to-pi+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, pi-, vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.vmsub.mults_piminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.pi-", comment="HERMES proton-to-pi- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=-1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, k+, vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.vmsub.mults_kplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.k+", comment="HERMES proton-to-k+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, k-, vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.vmsub.mults_kminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.k-", comment="HERMES proton-to-k- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=-2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, pi+, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.no-vmsub.mults_piplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.pi+.no-vmsub", comment="HERMES proton-to-pi+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, pi-, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.no-vmsub.mults_piminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.pi-.no-vmsub", comment="HERMES proton-to-pi- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=-1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, k+, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.no-vmsub.mults_kplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.k+.no-vmsub", comment="HERMES proton-to-k+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- proton target, k-, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.proton.zxpt-3D.no-vmsub.mults_kminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.p.k-.no-vmsub", comment="HERMES proton-to-k- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=1, h_2=-2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, pi+, vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.vmsub.mults_piplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.pi+", comment="HERMES deutron-to-pi+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, pi-, vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.vmsub.mults_piminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.pi-", comment="HERMES deutron-to-pi- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=-1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, k+, vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.vmsub.mults_kplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.k+", comment="HERMES deutron-to-k+ (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, k-, vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.vmsub.mults_kminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.k-", comment="HERMES deutron-to-k- (zxpt-3D) vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=-2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, pi+, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.no-vmsub.mults_piplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.pi+.no-vmsub", comment="HERMES deutron-to-pi+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, pi-, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.no-vmsub.mults_piminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.pi-.no-vmsub", comment="HERMES deutron-to-pi- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=-1, proc_id=proc_id, M_product=m_pion)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, k+, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.no-vmsub.mults_kplus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.k+.no-vmsub", comment="HERMES deutron-to-k+ (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- deuteron target, k-, no-vmsub
rows = _read_hermes_table(path_to_data + "hermes.deuteron.zxpt-3D.no-vmsub.mults_kminus.datnew_fortran")
ds = DataSet.empty("SIDIS", name="hermes3D.d.k-.no-vmsub", comment="HERMES deutron-to-k- (zxpt-3D) no-vmsub. thFactor contains DIS normalization (by MSHT20nnlo)",
                    reference="1212.5407", normErr=[], isNormalized=False)
_add_hermes_points(ds, rows, ps_def, h_1=12, h_2=-2, proc_id=proc_id, M_product=m_kaon)
ds.save_csv(path_to_save + ds.name + ".csv")
