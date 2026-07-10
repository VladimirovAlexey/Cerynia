"""
Parsing of E154 d2 moment data, neutron.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "E154" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- value hand-typed
directly from the paper text (no raw data file); pure reformat, nothing
changed (see RQCD.py / feedback_g2_conventions).

NOTE: the old script computes s using M_proton even though this is a
neutron target (`s=2*48.3*M_proton+M_proton**2`) -- ported verbatim, not
"fixed" to M_neutron, matching the old script's own (evidently intentional
or at least unchanged) convention.

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=101 (neutron).

NOTE: normErr=[].

Dataset (hep-ex/9705017, "text after (6)"), 1 point:
    E154_d2
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

s = 2 * 48.3 * M_proton + M_proton**2
Q_avg = sqrt(3.6)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

ds = DataSet.empty("D2", name="E154_d2",
                    comment="Taken from the text after (6)",
                    reference="hep-ex/9705017", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=-0.004, uncorrErr_0=0.038, uncorrErr_1=0.005))

ds.save_csv(path_to_save + ds.name + ".csv")
