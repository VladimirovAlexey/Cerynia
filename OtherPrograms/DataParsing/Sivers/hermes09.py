"""
Parsing of HERMES (2009 publication) transverse single-spin asymmetry
A_UT^Sivers for pi+/pi-/pi0/k+/k-, proton target: Q-integrated (dz/dx/dpt),
and Q<2/Q>2-sliced (pi+/k+ only, dz/dpt).

Source: /data/arTeMiDe_Repository/data/HERMES-SSA/ssa1.cgi (Q-integrated),
        ssa2.cgi (Q<2/Q>2)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/*.csv

Ported from the "hermes.sivers.*" (non-3d) blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadSSA_HERMES.py -- s, x/z/pT/Q
bins, cuts, and error treatment kept identical to the old parsing, except
one bugfix (see NOTE on Q>2 datasets below). The 3D (x,z,pT)-binned
datasets are a separate, newer measurement -- see hermes3D.py.

NOTE: process code. Confirmed correct as-is by the user -- primary process
[1,1,h_2,12001] and weight process [1,1,h_2,2001] (h_2: pi+=1, pi-=-1,
pi0=3, k+=2, k-=-2), matching the old script's "process"/"weightProcess"
fields exactly. h_1=1 (proton) is correct here (unlike the JLab case)
since HERMES genuinely uses a proton target.

NOTE: naming convention, given by the user. "hermes.sivers.*" ->
"hermes09.sivers.*" (2009 publication year, reference 0906.3918),
distinguishing these from the newer 3D data in hermes3D.py.

NOTE: Q>2 bugfix. The old script's Qbounds(xMin,xMax) helper (used for
every "Qint" dataset) takes a sqrt to convert Q^2 bounds into Q. But the 4
"Q>2" datasets (hermes09.sivers.pi+/k+.Q>2.dz/dpt) instead set
Q_max = xMax*yMax*sM2 directly, WITHOUT the sqrt -- giving a physically
implausible Q_max~19.7 GeV for a 27.6 GeV beam (vs ~4.4 GeV with the sqrt,
consistent with the <Q> values actually seen in this data, ~2-5 GeV).
Confirmed via direct comparison of raw <Q^2> values that "Q>2" is genuine,
independent data (not a duplicate/subset of "Qint" -- <Q^2>~5.5 for Q>2
vs ~2.4 for Qint), so this bug does matter. Per user confirmation, fixed:
Q_max = sqrt(xMax*yMax*sM2), matching the Qbounds() convention.

NOTE: x/z/pT bins are FIXED, shared bin-edge lists (xBin/zBin/ptBin below)
"presented by Gunar Schnell" -- for a "dz"/"dx"/"dpt" dataset, only the
named variable is picked per-point from its bin-edge list (by row/bin
number); the other two variables use their FULL range (bin-integrated).
Old dataset comment: "The data MUST be evaluated at a point" (i.e. the bin
integration is already folded into the systematic error) -- ported
verbatim including the comment, same pattern as jlab.sivers.* in JLab.py.

NOTE: ssa2.cgi rows have 2 extra trailing columns (VM-fraction and its
stat error) beyond the 9 columns shared with ssa1.cgi -- unused, ignored,
same as old parsing.

NOTE: thFactor=1 (ported verbatim; old code flags it "### tobe updated",
matches the established ratio/asymmetry-observable convention).

Datasets, 23 total:
    hermes09.sivers.{pi+,pi-,pi0,k+,k-}.Qint.{dz,dx,dpt}  -- 15 sets
    hermes09.sivers.{pi+,k+}.{Q<2,Q>2}.{dz,dpt}            -- 8 sets
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/HERMES-SSA/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/Sivers/"

M_proton = 0.938
m_pion = 0.139
m_kaon = 0.494

### 27.6 GeV beam on a fixed proton target (as in old parsing)
s = 2 * 27.6 * 0.938 + 0.938**2
### 7.3% overall scale uncertainty (as in old parsing)
normErr = [0.073]
includeCuts = True
### y, W^2 cuts (as in old parsing)
cutParams = [0.1, 0.95, 10., 10000.]

### bins presented by Gunar Schnell; the data MUST be evaluated at a point
### (bin integration is folded into the systematic error), see module NOTE
xBin = [0.023, 0.045, 0.067, 0.086, 0.113, 0.160, 0.220, 0.400]
zBin = [0.2, 0.27, 0.34, 0.41, 0.49, 0.56, 0.63, 0.70]
ptBin = [0.00001, 0.17, 0.25, 0.33, 0.41, 0.58, 0.80, 2.00]


def Qbounds(xMin, xMax):
    """
    Kinematic Q bounds given an x window (ported verbatim from old parsing):
    Qmin^2 = MAX(Q2min, xMin*yMin*(s-M^2), xMin/(1-xMin)*(W2min-M^2))
    Qmax^2 = MIN(Q2max, xMax*yMax*(s-M^2), xMax/(1-xMax)*(W2max-M^2))
    """
    Q2min, Q2max = 1., 10000.
    WM2min, WM2max = 10. - 0.938**2, 10000. - 0.938**2
    yMin, yMax = 0.1, 0.95
    sM2 = 2 * 27.6 * 0.938
    if xMax < 1:
        return [sqrt(max(Q2min, xMin * yMin * sM2, xMin / (1 - xMin) * WM2min)),
                sqrt(min(Q2max, xMax * yMax * sM2, xMax / (1 - xMax) * WM2max))]
    else:
        return [sqrt(max(Q2min, xMin * yMin * sM2, xMin / (1 - xMin) * WM2min)),
                sqrt(min(Q2max, xMax * yMax * sM2, 1000 * WM2max))]


def _read_hermes09_table(path, start, n_rows=7):
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    return [[float(v) for v in line.split()] for line in lines[start:start + n_rows]]


def _add_hermes09_points(ds, rows, diff_var, h_2, m_product, Q_min, Q_max_fn):
    """
    diff_var: which variable is binned per-row ("z", "x", or "pt"); the
    other two use their full range. Q_max_fn(x_bin) computes Q_max from the
    point's x-window (Qbounds() for Qint sets, or a fixed/derived value for
    the Q<2/Q>2 sets).
    """
    for i, row in enumerate(rows):
        binN = int(row[0])
        if diff_var == "z":
            x_bin, z_bin, pT_bin = (xBin[0], xBin[7]), tuple(zBin[binN - 1:binN + 1]), (ptBin[0], ptBin[7])
        elif diff_var == "x":
            x_bin, z_bin, pT_bin = tuple(xBin[binN - 1:binN + 1]), (zBin[0], zBin[7]), (ptBin[0], ptBin[7])
        else:  # "pt"
            x_bin, z_bin, pT_bin = (xBin[0], xBin[7]), (zBin[0], zBin[7]), tuple(ptBin[binN - 1:binN + 1])

        ds.add_point(dict(
            id=f"{ds.name}.{i}", ps_def=1, h_1=1, h_2=h_2, proc_id=12001,
            ps_def_weight=1, h_1_weight=1, h_2_weight=h_2, proc_id_weight=2001,
            s=s, Q_min=Q_min(x_bin) if callable(Q_min) else Q_min,
            Q_max=Q_max_fn(x_bin), Q_avg=sqrt(row[1]),
            x_min=x_bin[0], x_max=x_bin[1], x_avg=row[2],
            z_min=z_bin[0], z_max=z_bin[1], z_avg=row[4],
            pT_min=pT_bin[0], pT_max=pT_bin[1], pT_avg=row[5],
            M_target=M_proton, M_product=m_product,
            includeCuts=includeCuts,
            cutParams_0=cutParams[0], cutParams_1=cutParams[1],
            cutParams_2=cutParams[2], cutParams_3=cutParams[3],
            thFactor=1.,
            xSec=row[6],
            uncorrErr_0=row[7],
            uncorrErr_1=row[8],
        ))


def _Qmin_int(x_bin):
    return Qbounds(x_bin[0], x_bin[1])[0]


def _Qmax_int(x_bin):
    return Qbounds(x_bin[0], x_bin[1])[1]


### fixed vs old parsing (see module NOTE): Q_max = sqrt(xMax*yMax*sM2)
def _Qmax_hi(x_bin):
    return sqrt(2. * 27.6 * 0.938 * 0.95 * x_bin[1])


#%% -- pi+, Qint, dz
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 16)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Qint.dz",
                    comment="HERMES SSA-Sivers pi+ (integrated in Q, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 1, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, Qint, dx
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 28)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Qint.dx",
                    comment="HERMES SSA-Sivers pi+ (integrated in Q, differential in x). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "x", 1, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, Qint, dpt
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 40)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Qint.dpt",
                    comment="HERMES SSA-Sivers pi+ (integrated in Q, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 1, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi0, Qint, dz
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 52)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi0.Qint.dz",
                    comment="HERMES SSA-Sivers pi0 (integrated in Q, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 3, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi0, Qint, dx
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 64)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi0.Qint.dx",
                    comment="HERMES SSA-Sivers pi0 (integrated in Q, differential in x). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "x", 3, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi0, Qint, dpt
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 76)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi0.Qint.dpt",
                    comment="HERMES SSA-Sivers pi0 (integrated in Q, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 3, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, Qint, dz
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 88)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi-.Qint.dz",
                    comment="HERMES SSA-Sivers pi- (integrated in Q, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", -1, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, Qint, dx
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 100)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi-.Qint.dx",
                    comment="HERMES SSA-Sivers pi- (integrated in Q, differential in x). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "x", -1, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi-, Qint, dpt
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 112)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi-.Qint.dpt",
                    comment="HERMES SSA-Sivers pi- (integrated in Q, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", -1, m_pion, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Qint, dz
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 124)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Qint.dz",
                    comment="HERMES SSA-Sivers k+ (integrated in Q, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 2, m_kaon, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Qint, dx
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 136)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Qint.dx",
                    comment="HERMES SSA-Sivers k+ (integrated in Q, differential in x). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "x", 2, m_kaon, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Qint, dpt
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 148)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Qint.dpt",
                    comment="HERMES SSA-Sivers k+ (integrated in Q, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 2, m_kaon, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, Qint, dz
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 160)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k-.Qint.dz",
                    comment="HERMES SSA-Sivers k- (integrated in Q, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", -2, m_kaon, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, Qint, dx
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 172)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k-.Qint.dx",
                    comment="HERMES SSA-Sivers k- (integrated in Q, differential in x). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "x", -2, m_kaon, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k-, Qint, dpt
rows = _read_hermes09_table(path_to_data + "ssa1.cgi", 184)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k-.Qint.dpt",
                    comment="HERMES SSA-Sivers k- (integrated in Q, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", -2, m_kaon, _Qmin_int, _Qmax_int)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, Q<2, dz
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 16)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Q<2.dz",
                    comment="HERMES SSA-Sivers pi+ (Q<2, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 1, m_pion, 1., lambda xb: 2.)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, Q<2, dpt
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 28)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Q<2.dpt",
                    comment="HERMES SSA-Sivers pi+ (Q<2, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 1, m_pion, 1., lambda xb: 2.)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, Q>2, dz
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 40)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Q>2.dz",
                    comment="HERMES SSA-Sivers pi+ (Q>2, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 1, m_pion, 2., _Qmax_hi)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- pi+, Q>2, dpt
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 52)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.pi+.Q>2.dpt",
                    comment="HERMES SSA-Sivers pi+ (Q>2, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 1, m_pion, 2., _Qmax_hi)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Q<2, dz
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 64)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Q<2.dz",
                    comment="HERMES SSA-Sivers k+ (Q<2, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 2, m_kaon, 1., lambda xb: 2.)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Q<2, dpt
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 76)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Q<2.dpt",
                    comment="HERMES SSA-Sivers k+ (Q<2, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 2, m_kaon, 1., lambda xb: 2.)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Q>2, dz
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 88)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Q>2.dz",
                    comment="HERMES SSA-Sivers k+ (Q>2, differential in z). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "z", 2, m_kaon, 2., _Qmax_hi)
ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- k+, Q>2, dpt
rows = _read_hermes09_table(path_to_data + "ssa2.cgi", 100)
ds = DataSet.empty("SIDIS", name="hermes09.sivers.k+.Q>2.dpt",
                    comment="HERMES SSA-Sivers k+ (Q>2, differential in pt). The data MUST be evaluated at a point",
                    reference="0906.3918", normErr=normErr, isNormalized=False)
_add_hermes09_points(ds, rows, "pt", 2, m_kaon, 2., _Qmax_hi)
ds.save_csv(path_to_save + ds.name + ".csv")
