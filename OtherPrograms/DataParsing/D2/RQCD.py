"""
Parsing of RQCD (Regensburg) lattice-QCD d2 moment data.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "RQCD" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- all values are
hand-typed directly from the paper (no raw data file); pure reformat, the
old script's numbers are already correct/modern (per user instruction, same
exemption as G2 -- see feedback_g2_conventions), nothing changed here.

NOTE: process code. ps_def=1, h_1=1 fixed (per user instruction, same
convention as G2/D2 category); proc_id carried over unchanged from the old
single `process` integer. Unlike the rest of D2 (which reuses G2's
100/101/102=p/n/d scheme), these are lattice-QCD flavor-decomposition
results with no p/n/d target at all: RQCD_d2_ud uses proc_id=2 (point 0)
and proc_id=1 (point 1) -- the old script has NO inline comment saying
which of 1/2 is u vs d, so this is preserved exactly as coded, not
labelled. RQCD_d2_singlet uses proc_id=11 (u-d) and proc_id=12 (u+d), per
the old script's own inline comments. RQCD_d2_pn uses the ordinary
100=proton/101=neutron codes (a lattice p/n combination, "better not to use
because there are no sea-part" per the dataset's own comment).

NOTE: normErr=[] (matches every other D2/G2 dataset -- old normErr code was
always present but commented out).

Datasets (arXiv:2111.08306, table 4), 2 points each:
    RQCD_d2_ud, RQCD_d2_singlet, RQCD_d2_pn
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

ps_def, h_1 = 1, 1
Q_min, Q_max, Q_avg = 1.95, 2.05, 2.0
s = 4.0

#%% -- u and d quarks (no inline u/d label in the old script -- preserved as coded)
ds = DataSet.empty("D2", name="RQCD_d2_ud",
                    comment="Taken from table 4 for u and d quarks",
                    reference="arXiv:2111.08306", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=2,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=0.026, uncorrErr_0=0.004, uncorrErr_1=0.013))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=1,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=-0.0086, uncorrErr_0=0.0026, uncorrErr_1=0.0146))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- u-d and u+d quarks
ds = DataSet.empty("D2", name="RQCD_d2_singlet",
                    comment="Taken from table 4 for u-d and u+d quarks",
                    reference="arXiv:2111.08306", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=11,  # u-d
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=0.034, uncorrErr_0=0.004, uncorrErr_1=0.011))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=12,  # u+d
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=0.018, uncorrErr_0=0.005, uncorrErr_1=0.022))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- p and n combination
ds = DataSet.empty("D2", name="RQCD_d2_pn",
                    comment="Taken from table 4 for p and n combination (better not to use because there are no sea-part)",
                    reference="arXiv:2111.08306", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=0.0105, uncorrErr_0=0.019, uncorrErr_1=0.065))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=-0.0009, uncorrErr_0=0.0014, uncorrErr_1=0.0069))

ds.save_csv(path_to_save + ds.name + ".csv")
