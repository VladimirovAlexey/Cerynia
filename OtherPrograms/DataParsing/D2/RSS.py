"""
Parsing of RSS d2 moment data, proton/deuteron/neutron (two papers).

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "RSS (2006)" and "RSS (2010)" blocks of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- all values are
hand-typed directly from the paper text/tables (no raw data file); pure
reformat, nothing changed (see RQCD.py / feedback_g2_conventions).

NOTE: reference "hep-exp/0608003" (not the standard "hep-ex" prefix) and
"arxiv:0812:0031" (non-standard colon-separated arXiv id) are verbatim
typos in the old script, ported as-is, not corrected.

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=100/101/102
(proton/neutron/deuteron).

NOTE: both datasets share the same beam-energy-derived s
(`s=2*5.75*M_proton+M_proton**2`).

NOTE: RSS-2008_d2's 3rd point has uncorrErr_1=0.0 (a genuine
data-table value, not a normErr placeholder -- the "no zero normErr" rule
does not apply here since this is uncorrErr, not normErr).

NOTE: normErr=[].

Datasets:
    RSS-2006_d2 (hep-exp/0608003, "text in the last page"), 1 point (p)
    RSS-2008_d2 (arxiv:0812:0031, "table II", print statement calls this
        "RSS (2010)"), 3 points (p, d, n)
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

M_proton = 0.938
ps_def, h_1 = 1, 1

s = 2 * 5.75 * M_proton + M_proton**2

#%% -- RSS (2006): p
Q_avg = sqrt(1.3)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

ds = DataSet.empty("D2", name="RSS-2006_d2",
                    comment="Taken from text in the last page",
                    reference="hep-exp/0608003", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0057, uncorrErr_0=0.0009, uncorrErr_1=0.0007))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- RSS (2010) (table II): p, d, n
Q_avg = sqrt(1.28)
Q_min, Q_max = Q_avg - 0.5, Q_avg + 0.5

ds = DataSet.empty("D2", name="RSS-2008_d2",
                    comment="Taken from table II",
                    reference="arxiv:0812:0031", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0104, uncorrErr_0=0.0004, uncorrErr_1=0.0013))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=102,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0027, uncorrErr_0=0.0008, uncorrErr_1=0.0017))
ds.add_point(dict(id=f"{ds.name}.2", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=-0.0075, uncorrErr_0=0.0021, uncorrErr_1=0.0))

ds.save_csv(path_to_save + ds.name + ".csv")
