"""
Parsing of Hall A d2 moment data, neutron (three separate papers).

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "HallA-2004", "HallA-2014", "HallA-2016" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- all values are
hand-typed directly from the paper text/tables (no raw data file); pure
reformat, nothing changed (see RQCD.py / feedback_g2_conventions).

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=101 (neutron) for all
three datasets.

NOTE: all three share the same beam-energy-derived s
(`s=2*5.73*M_proton+M_proton**2`). HallA-2014_d2 and HallA-2016_d2 have
IDENTICAL xSec/Q values (HallA-2016 looks like an updated re-analysis of
the same measurement adding a 3rd systematic-error column) -- both kept as
separate datasets, per the old script (not merged/deduplicated).

NOTE: normErr=[].

Datasets:
    HallA-2004_d2 (hep-ex/0405006, "(24)"), 1 point
    HallA-2014_d2 (arxiv:1404.4003, "table I"), 2 points (2 Q-bins)
    HallA-2016_d2 (arxiv:1603.03612, "table X"), 2 points (2 Q-bins, 3rd
        uncorr. error column added vs. 2014)
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

s = 2 * 5.73 * M_proton + M_proton**2

#%% -- HallA-2004 ("(24)")
Q_avg = sqrt(5.0)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

ds = DataSet.empty("D2", name="HallA-2004_d2",
                    comment="Taken from (24)",
                    reference="hep-ex/0405006", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0062, uncorrErr_0=0.0028))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- HallA-2014 ("table I")
ds = DataSet.empty("D2", name="HallA-2014_d2",
                    comment="Taken from table I",
                    reference="arxiv:1404.4003", normErr=[], isNormalized=False)

Q_avg0 = sqrt(3.21)
Q_avg1 = sqrt(4.32)
ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_avg0 - 0.5, Q_max=Q_avg0 + 0.5, Q_avg=Q_avg0, thFactor=1.,
                   xSec=-0.00421, uncorrErr_0=0.00079, uncorrErr_1=0.00082))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_avg1 - 0.5, Q_max=Q_avg1 + 0.5, Q_avg=Q_avg1, thFactor=1.,
                   xSec=-0.00035, uncorrErr_0=0.00083, uncorrErr_1=0.00069))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- HallA-2016 ("table X")
ds = DataSet.empty("D2", name="HallA-2016_d2",
                    comment="Taken from table X",
                    reference="arxiv:1603.03612", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_avg0 - 0.5, Q_max=Q_avg0 + 0.5, Q_avg=Q_avg0, thFactor=1.,
                   xSec=-0.00421, uncorrErr_0=0.00079, uncorrErr_1=0.00082, uncorrErr_2=8e-05))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_avg1 - 0.5, Q_max=Q_avg1 + 0.5, Q_avg=Q_avg1, thFactor=1.,
                   xSec=-0.00035, uncorrErr_0=0.00083, uncorrErr_1=0.00069, uncorrErr_2=7e-05))

ds.save_csv(path_to_save + ds.name + ".csv")
