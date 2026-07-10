"""
Parsing of E155 d2 moment data, proton/deuteron (two separate papers).

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "E155-1999" and "E155" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- all values are
hand-typed directly from the paper text (no raw data file); pure reformat,
nothing changed (see RQCD.py / feedback_g2_conventions).

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=100/102 (proton/deuteron).

NOTE: both datasets use a fixed Q window of sqrt(5)+/-0.5 GeV, same pattern
as E143.py, but each has its own beam-energy-derived s.

NOTE: normErr=[].

Datasets, 2 points each (p, d):
    E155-1999_d2 (hep-ex/9901006, "text in page 5")
    E155_d2 (hep-ex/0204028, "text after (5)" -- print statement in the old
        script calls this "E155 (2002)")
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

Q_avg = sqrt(5.0)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

#%% -- E155-1999 (text in page 5): p, d
s = 2 * 38.8 * M_proton + M_proton**2

ds = DataSet.empty("D2", name="E155-1999_d2",
                    comment="Taken from the text in page 5",
                    reference="hep-ex/9901006", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.005, uncorrErr_0=0.008))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=102,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.008, uncorrErr_0=0.005))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E155 (2002) (text after (5)): p, d
s = 2 * 32.3 * M_proton + M_proton**2

ds = DataSet.empty("D2", name="E155_d2",
                    comment="Taken from the text after (5)",
                    reference="hep-ex/0204028", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0025, uncorrErr_0=0.0016, uncorrErr_1=0.001))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=102,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0054, uncorrErr_0=0.0023, uncorrErr_1=0.0005))

ds.save_csv(path_to_save + ds.name + ".csv")
