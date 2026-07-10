"""
Parsing of GHMP26 lattice-QCD d2 moment data (u-d combination).

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/D2/*.csv

Ported from the "GHMP26" block of
DataProcessor/OtherPrograms/DataParsingTw3/ParsingD2.py -- single value
hand-typed directly from the paper's equation (no raw data file); pure
reformat, nothing changed (see RQCD.py / feedback_g2_conventions).

NOTE: reference "arXiv:2604.00143" is a forward-dated arXiv identifier in
the old script -- ported verbatim, not "corrected".

NOTE: process code. ps_def=1, h_1=1 fixed; proc_id=11 (u-d, per the old
script's own inline comment, same code as RQCD_d2_singlet's u-d point --
but note thFactor differs: 1 here vs 2 there, preserved exactly as coded).

NOTE: normErr=[].

Dataset (arXiv:2604.00143, "eqn.(30).)" -- verbatim trailing typo preserved
from the old comment), 1 point:
    GHMP26_d2
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/D2/"

ps_def, h_1 = 1, 1
Q_min, Q_max, Q_avg = 1.95, 2.05, 2.0
s = 4.0

ds = DataSet.empty("D2", name="GHMP26_d2",
                    comment="Taken from eqn.(30).)",
                    reference="arXiv:2604.00143", normErr=[], isNormalized=False)

ds.add_point(dict(id=f"{ds.name}.0", ps_def=ps_def, h_1=h_1, proc_id=11,  # u-d
                   s=s, Q_min=Q_min, Q_max=Q_max, Q_avg=Q_avg, thFactor=1.,
                   xSec=0.0024, uncorrErr_0=0.0046))

ds.save_csv(path_to_save + ds.name + ".csv")
