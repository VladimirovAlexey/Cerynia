"""
Parsing of E615 pion-beam Drell-Yan (pi- + W -> mu+ mu- X), Q-differential and
xF-differential cross sections.

Source: /data/arTeMiDe_Repository/data/FNAL-615/FNAL-615(dQ).dat, FNAL-615(dx).dat
        (hand-typed tables, compiled from Conway et al. Phys.Rev.D 39 (1989))
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/piDY/*.csv

Ported from the "E615(pi)-dQ"/"E615(pi)-dxF" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadPionDYdataFiles.py -- s, Q/xF
ranges, error treatment, and the >0-uncertainty row filter are kept identical
to the old parsing. The old script's 20 per-single-bin child datasets
("E615(pi)-dQ-4.05", ..., "E615(pi)-dxF-0.9") are pure subsets of the full
dQ/dxF datasets below (same points, filtered to one bin) -- per user
instruction, these are NOT ported (redundant with the full datasets).

NOTE: process code [2,2,184074,1] -- same convention as E537.py (see that
file's NOTE for the h_2/proc_id modernization rationale).

NOTE: the DY Point schema's y_min/y_max fields hold the FEYNMAN-x (xF)
window here, not real rapidity -- same repurposing as E537.py. For "dQ",
y_min/y_max = [0.0, 1.0] is fixed (full integrated xF acceptance) for every
point; for "dxF", y_min/y_max is the actual per-point xF bin, and
Q_min/Q_max = [4.05, 8.55] is fixed instead.

NOTE: atmdeFactor. Same standing rule as E537.py, per the same user
confirmation -- the 3D DY atmdeFactor = (Q_max-Q_min)*(y_max-y_min)*
(qT_max-qT_min) is applied here too.

NOTE: raw qT is given as a bin-center value with a fixed +-0.125 GeV
half-width (not explicit low/high edges); the 0.001 factor in thFactor is a
pb->nb unit conversion (old comment), ported verbatim. Rows with a raw
uncertainty <=0 are dropped (as in old parsing).

Datasets (Phys.Rev.D 39 (1989) 92-122):
    E615.pi.dQ  -- 10 Q-bins (4.05-13.05 GeV), qT-differential, xF-integrated
    E615.pi.dxF -- 10 xF-bins (0.0-1.0), qT-differential, Q=[4.05,8.55] fixed
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/FNAL-615/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/piDY/"

### pi- beam on tungsten (W) target, modern nuclear code (per user instruction)
ps_def, h_1, h_2, proc_id = 2, 2, 184074, 1
### fixed 252 GeV/c pi- beam on a stationary nucleon target (as in old parsing)
s = 473.6
### 16% overall systematic (as in old parsing)
normErr = [0.16]
includeCuts = False


def _read_e615_table(path, tag_prefix):
    """
    Read a tab-separated E615 table, grouped into blocks by a "#<tag>
    [lo,hi]" header line. Returns a list of (bin_low, bin_high, rows),
    rows = float lists [pT_center, xSec, err].
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    blocks = []
    bin_range = None
    rows = None
    for line in lines:
        if line.startswith(tag_prefix):
            if rows is not None:
                blocks.append((bin_range[0], bin_range[1], rows))
            label = line[len(tag_prefix):].strip("[]")
            bin_range = [float(v) for v in label.split(",")]
            rows = []
        elif line == "" or line[0] == "#":
            continue
        else:
            rows.append([float(v) for v in line.split("\t")])
    if rows is not None:
        blocks.append((bin_range[0], bin_range[1], rows))
    return blocks


#%% -- Q-differential (E615(pi)-dQ)
blocks = _read_e615_table(path_to_data + "FNAL-615(dQ).dat", "#Q= ")

ds = DataSet.empty("DY", name="E615.pi.dQ", comment="E615 pi-data Q-differential",
                    reference="Phys.Rev.D 39 (1989) 92-122", normErr=normErr, isNormalized=False)

y_min, y_max = 0.0, 1.0  # full integrated xF acceptance, stashed in the y-fields
k = 0
for Q_min, Q_max, rows in blocks:
    for row in rows:
        qT_min, qT_max = row[0] - 0.125, row[0] + 0.125
        uncorrErr_0 = row[2]
        if uncorrErr_0 <= 0:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"{ds.name}.{k}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            ### pb->nb unit conversion (as in old parsing)
            thFactor=atmdeFactor / (qT_max - qT_min) / (Q_max - Q_min) * 0.001,
            includeCuts=includeCuts,
            xSec=row[1],
            uncorrErr_0=uncorrErr_0,
        ))
        k += 1

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- xF-differential (E615(pi)-dxF)
blocks = _read_e615_table(path_to_data + "FNAL-615(dx).dat", "#X= ")

ds = DataSet.empty("DY", name="E615.pi.dxF", comment="E615 pi-data x-differential",
                    reference="Phys.Rev.D 39 (1989) 92-122", normErr=normErr, isNormalized=False)

Q_min, Q_max = 4.05, 8.55  # fixed mass window
k = 0
for y_min, y_max, rows in blocks:
    for row in rows:
        qT_min, qT_max = row[0] - 0.125, row[0] + 0.125
        uncorrErr_0 = row[2]
        if uncorrErr_0 <= 0:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"{ds.name}.{k}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            ### pb->nb unit conversion (as in old parsing)
            thFactor=atmdeFactor / (qT_max - qT_min) / (y_max - y_min) * 0.001,
            includeCuts=includeCuts,
            xSec=row[1],
            uncorrErr_0=uncorrErr_0,
        ))
        k += 1

ds.save_csv(path_to_save + ds.name + ".csv")
