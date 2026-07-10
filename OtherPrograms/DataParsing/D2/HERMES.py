"""
Parsing of HERMES d2 moment data, proton.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "HERMES" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- value hand-typed
directly from the paper text (no raw data file); pure reformat, nothing
changed (see RQCD.py / feedback_g2_conventions).

NOTE: same underlying HERMES paper as G2's HERMES/HERMES.av datasets
(arXiv:1112.5584).

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=100 (proton).

NOTE: normErr=[].

Dataset (arxiv:1112.5584, "text in the last page"), 1 point:
    HERMES_d2
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

s = 2 * 27.6 * M_proton + M_proton**2
Q_avg = sqrt(5.0)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

ds = DataSet.empty("D2", name="HERMES_d2",
                    comment="Taken from text in the last page",
                    reference="arxiv:1112.5584", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0148, uncorrErr_0=0.0096, uncorrErr_1=0.0048))

ds.save_csv(path_to_save + ds.name + ".csv")
