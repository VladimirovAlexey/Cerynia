"""
Parsing of JLab (Hall A, polarized 3He / effective-neutron target) transverse
single-spin asymmetry A_UT^Sivers for pi+/pi-/k+/k-.

Source: /data/arTeMiDe_Repository/data/JLab-Sivers/Siv_{pi+,pi-,k+,k-}
        (hand-typed tables; comma decimal separator, tab-separated columns)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from the "jlab.sivers.*" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSSA_JLAB.py -- s, x/z/pT/Q bins,
cuts, and error treatment kept identical to the old parsing.

NOTE: process code. Old parsing used [1,1,h_2,12003] (primary) /
[1,1,h_2,2003] (weight), commented "neutron target" despite h_1=1 (which in
the modern SIDIS convention, see HERMES3D.py/COMPASS.py, means a PROTON
target -- h_1=12 is deuteron). The raw data files themselves also explicitly
label "target: neutron". Per user instruction, updated to h_1=11 (neutron)
and primary proc_id=12001 (replacing 12003); weight proc_id_weight=2001
(the standing SIDIS Sivers weight code, same as every other SIDIS sub-case
in this category, replacing the old scheme's per-hadron 2003/2001/2002...).
h_2 unchanged: 1=pi+, -1=pi-, 2=k+, -2=k-.

NOTE: M_target is still the proton mass (0.938), ported verbatim -- old
parsing already used the proton mass even under its own "neutron target"
label (proton/neutron mass difference is ~0.06%, not corrected here, matches
established "port verbatim" policy).

NOTE: x/z/pT/Q bins are FIXED, shared, hand-eyeballed ranges from the papers
("bins are looked by eye from the 1404.7204 paper"), not real per-point bin
edges -- the old dataset comment explicitly flags "The data MUST be
evaluated at a point" (i.e. artemide should evaluate theory at the <x>/<z>/
<pT>/<Q> central values, not integrate over these bins). Ported verbatim,
including the comment.

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated", but
matches the established ratio/asymmetry-observable convention from
DY_angular/Sivers-STAR -- artemide builds the weighted asymmetry via the
weight-process columns).

NOTE: raw files use a comma decimal separator and inconsistent tab counts
(only Siv_pi+ has trailing garbage tabs); a shared reader normalizes both.

Datasets:
    jlab.sivers.pi+ -- 4 points, reference 1106.0363
    jlab.sivers.pi- -- 4 points, reference 1106.0363
    jlab.sivers.k+  -- 4 points, reference 1404.7204
    jlab.sivers.k-  -- 1 point,  reference 1404.7204
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/JLab-Sivers/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_proton = 0.938
m_pion = 0.139
m_kaon = 0.494

### 5.9 GeV electron beam on a fixed nucleon target (as in old parsing)
s = 2 * 5.9 * 0.938 + 0.938**2
### hand-eyeballed shared bins (not real per-point edges, see module NOTE)
x_min, x_max = 0.1, 0.4
z_min, z_max = 0.45, 0.6
pT_min, pT_max = 0.1, 0.6
Q_min, Q_max = 1., 1.7
includeCuts = False
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.95, 2.3, 10000.]

### neutron target (per user instruction and raw files' own "target: neutron" label)
ps_def, h_1 = 1, 11
proc_id_weight = 2001


def _read_jlab_table(path, n_rows):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    rows = []
    for line in lines[1:1 + n_rows]:
        line = line.replace(",", ".").replace("\t\t\t", "")
        rows.append([float(v) for v in line.split("\t")])
    return rows


#%% -- pi+
rows = _read_jlab_table(path_to_data + "Siv_pi+", 4)

ds = DataSet.empty("SIDIS", name="jlab.sivers.pi+",
                    comment="JLab SSA-Sivers pi+. The data MUST be evaluated at a point",
                    reference="1106.0363", normErr=[], isNormalized=False)

for i, row in enumerate(rows):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=1, proc_id=12001,
        ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=1, proc_id_weight=proc_id_weight,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[4]),
        x_min=x_min, x_max=x_max, x_avg=row[1],
        z_min=z_min, z_max=z_max, z_avg=row[3],
        pT_min=pT_min, pT_max=pT_max, pT_avg=row[5],
        M_target=M_proton, M_product=m_pion,
        thFactor=1.,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=row[6],
        uncorrErr_0=row[7],
        uncorrErr_1=row[8],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-
rows = _read_jlab_table(path_to_data + "Siv_pi-", 4)

ds = DataSet.empty("SIDIS", name="jlab.sivers.pi-",
                    comment="JLab SSA-Sivers pi-. The data MUST be evaluated at a point",
                    reference="1106.0363", normErr=[], isNormalized=False)

for i, row in enumerate(rows):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=-1, proc_id=12001,
        ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=-1, proc_id_weight=proc_id_weight,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[4]),
        x_min=x_min, x_max=x_max, x_avg=row[1],
        z_min=z_min, z_max=z_max, z_avg=row[3],
        pT_min=pT_min, pT_max=pT_max, pT_avg=row[5],
        M_target=M_proton, M_product=m_pion,
        thFactor=1.,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=row[6],
        uncorrErr_0=row[7],
        uncorrErr_1=row[8],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+
rows = _read_jlab_table(path_to_data + "Siv_k+", 4)

ds = DataSet.empty("SIDIS", name="jlab.sivers.k+",
                    comment="JLab SSA-Sivers k+. The data MUST be evaluated at a point",
                    reference="1404.7204", normErr=[], isNormalized=False)

for i, row in enumerate(rows):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=2, proc_id=12001,
        ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=2, proc_id_weight=proc_id_weight,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[4]),
        x_min=x_min, x_max=x_max, x_avg=row[1],
        z_min=z_min, z_max=z_max, z_avg=row[3],
        pT_min=pT_min, pT_max=pT_max, pT_avg=row[5],
        M_target=M_proton, M_product=m_kaon,
        thFactor=1.,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        ### columns 6,7 are W,W' (not used); value/stat/systabs are 8,9,11 (systrel at 10 unused)
        xSec=row[8],
        uncorrErr_0=row[9],
        uncorrErr_1=row[11],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-
rows = _read_jlab_table(path_to_data + "Siv_k-", 1)

ds = DataSet.empty("SIDIS", name="jlab.sivers.k-",
                    comment="JLab SSA-Sivers k-. The data MUST be evaluated at a point",
                    reference="1404.7204", normErr=[], isNormalized=False)

for i, row in enumerate(rows):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=ps_def, h_1=h_1, h_2=-2, proc_id=12001,
        ps_def_weight=ps_def, h_1_weight=h_1, h_2_weight=-2, proc_id_weight=proc_id_weight,
        s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=sqrt(row[4]),
        x_min=x_min, x_max=x_max, x_avg=row[1],
        z_min=z_min, z_max=z_max, z_avg=row[3],
        pT_min=pT_min, pT_max=pT_max, pT_avg=row[5],
        M_target=M_proton, M_product=m_kaon,
        thFactor=1.,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=row[8],
        uncorrErr_0=row[9],
        uncorrErr_1=row[11],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
