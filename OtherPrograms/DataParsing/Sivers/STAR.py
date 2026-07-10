"""
Parsing of STAR transverse single-spin asymmetry A_N for W+/W- (differential
in rapidity) and Z (single point) boson production.

Source: hardcoded literal values in the old parsing script (see NOTE below) --
        no raw data file is read for these 3 datasets.
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from the "star26.sivers.W+.dy"/"star26.sivers.W-.dy"/"star23.sivers.Z"
blocks of the old DataProcessor/OtherPrograms/DataParsing/ReadSSA_STAR_update.py
-- s, Q/y/qT ranges, cuts, and error treatment kept identical to the old
parsing, except one bugfix (see NOTE on star26.sivers.W-.dy below).

NOTE: process code. Confirmed correct as-is by the user -- primary process
[1,1,1,10003/10004/10005] (Z/W+/W-) and weight process [1,1,1,3/4/5],
matching the old script's "process"/"weightProcess" fields exactly. The
weight-process columns (ps_def_weight, h_1_weight, h_2_weight,
proc_id_weight) are the same mechanism introduced for DY_angular ratio
observables -- A_N is itself a ratio (spin asymmetry), so it needs the same
denominator-process bookkeeping.

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
but this matches the established ratio-observable convention from
DY_angular -- artemide builds the correctly-weighted asymmetry via the
weight-process columns, no separate Jacobian needed).

NOTE: star26.sivers.W+.dy and W-.dy have no real reference yet (old code:
reference="????.????", i.e. unpublished/private-communication data) --
ported verbatim.

Datasets:
    star23.sivers.Z    -- 1 point,  reference 2308.15496
    star26.sivers.W+.dy -- 3 points, reference unpublished (private communication)
    star26.sivers.W-.dy -- 3 points, reference unpublished (private communication)
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_Z = 91.19
M_W = 80.38

### beam-polarization scale uncertainty, common to all STAR A_N sets (as in old parsing)
normErr = [0.014]
s = 500.**2
includeCuts = False
### pT1, pT2, yMin, yMax (as in old parsing)
cutParams = [25., 25., -1., 1.]

#%% -- Z, single point
ds = DataSet.empty("DY", name="star23.sivers.Z", comment="STAR AN for Z",
                    reference="2308.15496", normErr=normErr, isNormalized=False)

ds.add_point(dict(
    id="star23.sivers.Z.0", ps_def=1, h_1=1, h_2=1, proc_id=10003,
    ps_def_weight=1, h_1_weight=1, h_2_weight=1, proc_id_weight=3,
    s=s, Q_min=73., Q_max=114., Q_avg=M_Z, y_min=-1., y_max=1.,
    qT_min=0.5, qT_max=10.,
    thFactor=1.,
    includeCuts=includeCuts,
    cutParams_0=cutParams[0], cutParams_1=cutParams[1],
    cutParams_2=cutParams[2], cutParams_3=cutParams[3],
    xSec=0.056,
    uncorrErr_0=0.081,
    uncorrErr_1=0.050,
))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- W+, differential in y (data sent by Oleg Eyser, private communication)
### original table, as received and copied verbatim into the old parsing script:
### W-plus
### y_reco   AN      stat    syst
### -0.50   -0.0539  0.0672  0.0289
###  0.00    0.0713  0.0474  0.0435
###  0.50   -0.1094  0.0670  0.0319
xSec_Wp = [-0.0539, 0.0713, -0.1094]
stat_Wp = [0.0672, 0.0474, 0.0670]
syst_Wp = [0.0289, 0.0435, 0.0319]
y_bins_Wp = [[-0.7, -0.3], [-0.3, 0.3], [0.3, 0.7]]

ds = DataSet.empty("DY", name="star26.sivers.W+.dy", comment="STAR AN for W+ (differential in y)",
                    reference="????.????", normErr=normErr, isNormalized=False)

for i, (y_min, y_max) in enumerate(y_bins_Wp):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=1, h_1=1, h_2=1, proc_id=10004,
        ps_def_weight=1, h_1_weight=1, h_2_weight=1, proc_id_weight=4,
        s=s, Q_min=M_W - 30., Q_max=M_W + 30., Q_avg=M_W, y_min=y_min, y_max=y_max,
        qT_min=0.5, qT_max=10.,
        thFactor=1.,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec_Wp[i],
        uncorrErr_0=stat_Wp[i],
        uncorrErr_1=syst_Wp[i],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- W-, differential in y (data sent by Oleg Eyser, private communication)
### original table, as received and copied verbatim into the old parsing script:
### W-minus
### y_reco   AN      stat    syst
### -0.47   -0.0618  0.1426  0.0692
### -0.01   -0.0212  0.0965  0.0302
###  0.47    0.0712  0.1517  0.0242
xSec_Wm = [-0.0618, -0.0212, 0.0712]
stat_Wm = [0.1426, 0.0965, 0.1517]
syst_Wm = [0.0692, 0.0302, 0.0242]
y_bins_Wm = [[-0.7, -0.3], [-0.3, 0.3], [0.3, 0.7]]

ds = DataSet.empty("DY", name="star26.sivers.W-.dy", comment="STAR AN for W- (differential in y)",
                    reference="????.????", normErr=normErr, isNormalized=False)

for i, (y_min, y_max) in enumerate(y_bins_Wm):
    ds.add_point(dict(
        id=f"{ds.name}.{i}", ps_def=1, h_1=1, h_2=1, proc_id=10005,
        ps_def_weight=1, h_1_weight=1, h_2_weight=1, proc_id_weight=5,
        s=s, Q_min=M_W - 30., Q_max=M_W + 30., Q_avg=M_W, y_min=y_min, y_max=y_max,
        qT_min=0.5, qT_max=10.,
        thFactor=1.,
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec_Wm[i],
        uncorrErr_0=stat_Wm[i],
        uncorrErr_1=syst_Wm[i],
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
