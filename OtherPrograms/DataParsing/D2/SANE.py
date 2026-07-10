"""
Parsing of SANE d2 moment data, proton.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "SANE" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- values hand-typed
directly from the paper's table 1 (no raw data file); pure reformat,
nothing changed (see RQCD.py / feedback_g2_conventions).

NOTE: unlike every other D2 dataset, each point has its OWN Q window and
beam energy (no shared fixed-window/fixed-s pattern): point 0 uses a 4.7
GeV beam and Q=[sqrt(2),sqrt(3.5)]; point 1 uses a 5.9 GeV beam and
Q=[sqrt(3.5),sqrt(5)] -- both Q windows are "guessed by us" per the old
comment (no HEPData Q-range given by the experiment).

NOTE: thFactor=-1 (unique sign flip in this category -- the paper defines
d2 with an extra factor of -1 vs. the standard convention, see their
eqn.(2), per the old comment).

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=100 (proton).

NOTE: normErr=[].

Dataset (arXiv:1805.08835, "table 1"), 2 points:
    SANE_d2
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

ds = DataSet.empty("D2", name="SANE_d2",
                    comment="Taken from table 1 (they have factor -1 in definition, see eqn.(2); Q-range is guessed by us)",
                    reference="arXiv:1805.08835", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=2 * M_proton * 4.7 + M_proton**2,
                   Q_min=sqrt(2.), Q_max=sqrt(3.5), Q_avg=sqrt(2.8), thFactor=-1.,
                   xSec=-0.00414, uncorrErr_0=0.00205, uncorrErr_1=0.00256))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=2 * M_proton * 5.9 + M_proton**2,
                   Q_min=sqrt(3.5), Q_max=sqrt(5.), Q_avg=sqrt(4.3), thFactor=-1.,
                   xSec=-0.00149, uncorrErr_0=0.00156, uncorrErr_1=0.00368))

ds.save_csv(path_to_save + ds.name + ".csv")
