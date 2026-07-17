"""
Demonstrations of some elemental operations with Cerynia
In this example, it is shown
1) how to create a DataSet from scratch (no reading from a csv file)
2) which fields of a point are mandatory, and which are optional
3) what is thFactor, and how normalization/normalization-uncertainty works
4) how to save the resulting DataSet to a csv file (in the DataLib format)
"""
#######################################
# The first part of any code is to load libraries
# In this example we need only Cerynia (no artemide/harpy is used here,
# since we only create and save data, we do not compute any theory)
#######################################

import os
import Cerynia

### this is just the current folder; we will save the example csv-file here
THIS_DIR = os.path.dirname(__file__)

#%%
### A DataSet always belongs to a certain processType (see Cerynia.Point.PROCESS_TYPES).
### Currently defined types are "DY" (Drell-Yan), "SIDIS", "G2" and "D2".
### The processType fixes which columns (fields) a point must have, because different
### processes are described by different kinematic variables.
###
### You can always print the schema (mandatory + optional columns) of a processType:
Cerynia.Point.schema("DY")
### try also Cerynia.Point.schema("SIDIS") to see the (different) SIDIS fields.

#%%
### Here we create a small DY DataSet, with point values inspired by the real
### CDF1 measurement (DataLib/DY/CDF1.csv), but hand-typed instead of read from file,
### to keep the example as simple as possible.
###
### DataSet.empty(...) creates a DataSet with no points yet; points are added one-by-one
### afterward via .add_point(). The arguments of empty() are the SET-level (not per-point)
### properties:
###   processType         : "DY" here (mandatory)
###   name, comment, reference : free-text bookkeeping (all optional, default "")
###   normErr              : list of GLOBAL normalization uncertainties, given as
###                          RELATIVE values (e.g. 0.039 = 3.9%). These are uncertainties
###                          common to ALL points of the set (e.g. luminosity uncertainty),
###                          as opposed to per-point uncorrErr_i/corrErr_i (see below).
###                          Optional, default: [] (no normalization uncertainty).
###   isNormalized         : True if the data is a NORMALIZED (shape-only) measurement,
###                          i.e. the overall scale of the measurement is not meaningful,
###                          and theory must be rescaled to match the data before chi2.
###                          Optional, default: False (CDF1 is an absolute cross-section
###                          measurement, so we set it to False here).
###   normalizationMethod  : how the rescaling of isNormalized=True data is done
###                          ("integral" or "bestChi2"), only relevant if isNormalized=True.
###                          Optional, default: "integral".
ds = Cerynia.DataSet.empty(
    "DY",
    name="CDF1-example",
    comment="A minimal hand-made subset of CDF1, for demonstration only",
    reference="hep-ex/0001021",
    normErr=[0.039],
    isNormalized=False,
)

#%%
### Now we populate the set with points. Each point is just a dict of column->value,
### passed to ds.add_point(point). Internally it is validated against the "DY" schema
### (Cerynia.Point.validate), which:
###   1) checks that all MANDATORY columns are present (raises ValueError if not)
###   2) fills all missing OPTIONAL columns with their defaults
###   3) checks that bin orderings are sane (Q_min<=Q_max, qT_min<=qT_max, etc.)
###
### MANDATORY fields for "DY" (see Cerynia.Point.REQUIRED["DY"]):
###   id        : a unique string identifying the point within the set
###   xSec      : the measured value (here, a cross-section)
###   s         : Mandelstam s = collision energy^2 (GeV^2)
###   Q_min/Q_max            : the invariant-mass bin of the lepton pair
###   y_min/y_max             : the rapidity bin of the lepton pair
###   qT_min/qT_max           : the transverse-momentum bin of the lepton pair
###   ps_def, h_1, h_2, proc_id : integers that tell artemide WHICH process/hadrons this is
###                              (e.g. h_1=1,h_2=-1 for proton-antiproton, as in CDF1;
###                              see Cerynia/harpyInterface.py and artemide-manual for the
###                              full list of codes -- these are NOT something you should
###                              invent yourself, always take them from a similar case)
###   thFactor  : a multiplicative theory factor, see explanation below
###   includeCuts : True/False, whether artemide should apply an internal kinematic cut
###                 to the point (defined by cutParams_0..3, see below) -- this is a cut
###                 applied INSIDE artemide's computation, unrelated to DataSet.cut()
###                 (which is a Cerynia-side filter/selection of points).
###
### OPTIONAL fields for "DY" (filled automatically if omitted):
###   Q_avg, y_avg, qT_avg : bin centers, computed as (min+max)/2 if not given
###   cutParams_0..3       : parameters of the includeCuts cut window;
###                          default (0.0, 0.0, -100.0, 100.0) is a very wide, effectively
###                          "no cut" window -- only meaningful when includeCuts=True
###   uncorrErr_0, uncorrErr_1, ...   : per-point UNCORRELATED errors (one random shift per
###                                    point). You can add as many uncorrErr_i columns as
###                                    you need; if a point/set has none, just omit them.
###   corrErr_0, corrErr_1, ...       : per-point CORRELATED errors (one shared random shift
###                                    per source, across all points of the set). Also fully
###                                    optional, add only if the measurement has them.
###
### NOTE on thFactor: artemide computes a bin-AVERAGED cross-section over (Q,y,qT).
### Experiments, however, may report something else (e.g. a bin-INTEGRATED cross-section,
### or use a different unit convention). thFactor is the number that artemide's raw
### result must be multiplied by to match what the experiment reports; the result of the
### computation is (theory from artemide) * thFactor, matched against xSec.
### For CDF1, the data is bin-integrated in (Q,y,qT), so thFactor equals the bin-volume
### (Q_max-Q_min)*(y_max-y_min)*(qT_max-qT_min) times an overall unit/convention factor
### fixed by comparison with the original publication (here, thFactor=1000. for all points,
### since all three example points share the same bin widths).
ds.add_point(dict(
    id="CDF1-example.0", ps_def=1, h_1=1, h_2=-1, proc_id=3,
    s=3240000.0, Q_min=66.0, Q_max=116.0, y_min=-10.0, y_max=10.0,
    qT_min=0.0, qT_max=0.5, thFactor=1000.0, includeCuts=False,
    xSec=3.35, uncorrErr_0=0.54,
))
ds.add_point(dict(
    id="CDF1-example.1", ps_def=1, h_1=1, h_2=-1, proc_id=3,
    s=3240000.0, Q_min=66.0, Q_max=116.0, y_min=-10.0, y_max=10.0,
    qT_min=0.5, qT_max=1.0, thFactor=1000.0, includeCuts=False,
    xSec=10.1, uncorrErr_0=1.0,
))
ds.add_point(dict(
    id="CDF1-example.2", ps_def=1, h_1=1, h_2=-1, proc_id=3,
    s=3240000.0, Q_min=66.0, Q_max=116.0, y_min=-10.0, y_max=10.0,
    qT_min=1.0, qT_max=1.5, thFactor=1000.0, includeCuts=False,
    xSec=14.8, uncorrErr_0=1.2,
))

### At this point the set already has all 3 points; we can inspect it
print(ds)
ds.info()

#%%
### RECUP: to build a DataSet by hand
### 1) Cerynia.DataSet.empty(processType, name=..., normErr=[...], isNormalized=...)
### 2) ds.add_point({...}) for each point (mandatory fields must be given;
###    optional fields fall back to their defaults; see Cerynia.Point.schema(processType))

#%%
### Finally, the DataSet can be saved to a csv file, in the same format used by DataLib.
### Set-level metadata (name, comment, reference, processType, isNormalized,
### normalizationMethod, normErr) is written as leading "# key: value" comment lines,
### followed by the points themselves as a normal csv table.
OUT_PATH = THIS_DIR + "/CreateDataFile_example.csv"
ds.save_csv(OUT_PATH)
print("Saved example DataSet to:", OUT_PATH)

### Now, this .csv is ready for operation!