"""
Parsing of QCDSF lattice-QCD d2 moment data, proton and neutron.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "QCDSF" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- both values are
hand-typed directly from the paper text (no raw data file); pure reformat,
nothing changed (see RQCD.py / feedback_g2_conventions for the category-wide
exemption from the audit policy).

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=100/101 (proton/neutron,
same convention as G2's target-species code).

NOTE: normErr=[].

Dataset (arXiv:2408.03621, "text in page 3.)" -- verbatim trailing typo
preserved from the old comment), 2 points:
    QCDSF_d2
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

ps_def, h_1 = 1, 1
Q_min, Q_max, Q_avg = 1.95, 2.05, 2.0
s = 4.0

ds = DataSet.empty("D2", name="QCDSF_d2",
                    comment="Taken from text in page 3.)",
                    reference="arXiv:2408.03621", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=100,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=0.046, uncorrErr_0=0.007, uncorrErr_1=0.016))
ds.add_point(dict(id=f"{ds.name}.1", ps_def=ps_def, h_1=h_1, proc_id=101,
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=2.,
                   xSec=0.023, uncorrErr_0=0.005, uncorrErr_1=0.008))

ds.save_csv(path_to_save + ds.name + ".csv")
