"""
Parsing of STAR 510 GeV (p+p) Z/gamma* -> l+l- transverse-momentum spectrum.

Source: /data/arTeMiDe_Repository/data/STAR-DY/HEPData-ins2692202-v1-Fig2.csv
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/STAR510.csv

Updated to read the official HEPData release (ins2692202, Fig2) instead of the
earlier hand-relayed STAR_data.txt -- the user confirmed this file is correct
and supersedes it. Only qT_avg, xSec, and the errors change; everything else
(process code, s, Q/y-range, cuts, normErr, the hardcoded qT-bin edges) is
kept identical to the old parsing / previous version of this script, cross-
checked against STAR-DY/comment-to-data.txt (fiducial cuts pT_lepton>25 GeV,
|eta_lepton|<1, 73<M<114 GeV, 8.5% luminosity -- all already matched).

NOTE: this file only gives each point's measured pT centroid, not the bin
edges -- reusing the same hardcoded bin-edge list as before
(bins=[0.,1.25,2.5,...]), ported verbatim from old parsing; not re-derivable
from the data alone. The measured centroid is kept as a qT_avg override
rather than defaulting to the bin midpoint, matching old parsing's explicit
p["<qT>"]=... assignment.

Raw format: HEPData CSV, 4 "#:" header lines + 1 column-header line, 11 data
rows, trailing blank line. Columns [p_T, differential cross section, stat +,
stat -, sys +, sys -] -- signed +/- pairs, symmetrized the usual way.

Dataset:
    STAR510 -- 2308.15496
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

M_Z = 91.  # forces <Q> = Z-boson mass rather than the true bin midpoint (93.5),
           # matching the old parsing's explicit p["<Q>"]=M_Z override

### qT bin edges (as in old parsing) -- the raw table only gives a measured
### centroid per point, not the bin edges themselves
_PT_BIN_EDGES = [0., 1.25, 2.5, 3.75, 5., 7.5, 10., 12.5, 15., 17.5, 20., 25.]

#%%
# ============================================================================
# STAR 510 GeV -- 2308.15496 (HEPData ins2692202, Fig2)
# ============================================================================
with open(path_to_data + "STAR-DY/HEPData-ins2692202-v1-Fig2.csv") as f:
    lines = [line.rstrip("\n") for line in f]
lines = lines[5:-1]  # drop 4 "#:" lines + column header, and the trailing blank line
rows = [[float(v) for v in line.split(",")] for line in lines]

Q_min, Q_max = 73., 114.
### Process is Z-DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 3
### full rapidity window (as in old parsing; no symmetrize factor needed)
y_min, y_max = -1., 1.
includeCuts = True
cutParams = [25., 25., -1., 1.]
### 8.5% normalization uncertainty (luminosity, per comment-to-data.txt; as in old parsing)
normErr = [0.085]

ds = DataSet.empty("DY", name="STAR510", comment="STAR DY data", reference="2308.15496",
                    normErr=normErr, isNormalized=False)

for i, (pt_centroid, xSec, statp, statm, sysp, sysm) in enumerate(rows):
    qT_min, qT_max = _PT_BIN_EDGES[i], _PT_BIN_EDGES[i + 1]

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"STAR510.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=510.**2, Q_min=Q_min, Q_max=Q_max, Q_avg=M_Z, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max, qT_avg=pt_centroid,
        ### divide by bin size, no symmetrize factor (as in old parsing)
        thFactor=atmdeFactor / (qT_max - qT_min),
        includeCuts=includeCuts,
        cutParams_0=cutParams[0], cutParams_1=cutParams[1],
        cutParams_2=cutParams[2], cutParams_3=cutParams[3],
        xSec=xSec,
        uncorrErr_0=(statp - statm) / 2.,
        uncorrErr_1=(sysp - sysm) / 2.,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
