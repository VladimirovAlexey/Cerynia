"""
Parsing of E143 d2 moment data, proton/deuteron/neutron.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "E143-1995" and "E143" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- all values are
hand-typed directly from the paper text/tables (no raw data file); pure
reformat, nothing changed (see RQCD.py / feedback_g2_conventions).

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=100/101/102
(proton/neutron/deuteron, same convention as G2).

NOTE: both datasets share the same beam-energy-derived s and a fixed
+/-0.5 Q window around sqrt(5), ported verbatim from the old script's
`s_current=2*29.*M_proton+M_proton**2`, `<Q>=sqrt(5.0)`,
`Q=[<Q>-0.5,<Q>+0.5]` pattern.

NOTE: normErr=[].

Datasets, 1 point each:
    E143-1995_d2 (hep-ex/9511013, "table 5a"), 2 points (p, d)
    E143_d2 (hep-ph/9802357, "table XXXIII" -- same underlying paper as
        G2's E143.p/d/n), 3 points (p, d, n)
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

s = 2 * 29. * M_proton + M_proton**2
Q_avg = sqrt(5.0)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

#%% -- E143-1995 (table 5a): p, d
ds = DataSet.empty("D2", name="E143-1995_d2",
                    comment="Taken from table 5a",
                    reference="hep-ex/9511013", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0054, uncorrErr_0=0.005))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=102,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0039, uncorrErr_0=0.0092))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- E143 (table XXXIII): p, d, n
ds = DataSet.empty("D2", name="E143_d2",
                    comment="Taken from table XXXIII",
                    reference="hep-ph/9802357", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0058, uncorrErr_0=0.005))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=102,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0051, uncorrErr_0=0.0092))
ds.add_point(dict(id=f"{ds.name}.2", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0050, uncorrErr_0=0.021))

ds.save_csv(path_to_save + ds.name + ".csv")
