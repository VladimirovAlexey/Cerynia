"""
Parsing of E537 pion-beam Drell-Yan (pi- + W -> mu+ mu- X), Q-differential and
xF-differential cross sections.

Source: /data/arTeMiDe_Repository/data/FNAL-537/piminus+W(dQ).dat, pminus+W(dxF).dat
        (HEPData-style hand-typed tables, hepdata.23243.v1 Tables 5 and 6)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/piDY/*.csv

Ported from the "E537(pi)-dQ"/"E537(pi)-dxF" blocks of the old
DataProcessor/OtherPrograms/DataParsing/ReadPionDYdataFiles.py -- s, Q/xF
ranges, error treatment, and the >0-uncertainty row filter are kept identical
to the old parsing.

NOTE: process code. Old parsing used [2,2,1,103] (h_2=1, not a real target
code). Per user instruction, updated to [2,2,184074,1]: h_2=184074 is the
modern nuclear code for the tungsten (W) target (matching the regular DY
branch's E537 sets); proc_id=1 is the modern generic/continuum DY code
(replacing the old scheme's 103).

NOTE: the DY Point schema's y_min/y_max fields hold the FEYNMAN-x (xF)
window here, not real rapidity -- ported verbatim from old parsing (which
did the same repurposing, no dedicated xF field exists). For "dQ", y_min/
y_max = [-0.1, 1.0] is fixed (the full integrated xF acceptance) for every
point; for "dxF", y_min/y_max is the actual per-point xF bin from the raw
table, and Q_min/Q_max = [4.0, 9.0] is fixed instead.

NOTE: atmdeFactor. Per user confirmation, the standard 3D DY atmdeFactor =
(Q_max-Q_min)*(y_max-y_min)*(qT_max-qT_min) is applied here too (multiplying
the old thFactor expression), even though the y-dimension is xF rather than
rapidity: "It is the definition (variable y for dxF experiments keeps xF)".

NOTE: raw values are pT^2 (not pT); qT_min/qT_max = sqrt(pT^2_low/high). The
0.001 factor in thFactor is a pb->nb unit conversion (old comment), ported
verbatim. Rows with a symmetrized uncertainty <=0 are dropped (as in old
parsing) -- a data-validity filter, not a physics redefinition.

Datasets (Phys.Rev.D 93 (1988) 1377):
    E537.pi.dQ  -- 10 Q-bins (4.0-9.0 GeV), qT-differential, xF-integrated
    E537.pi.dxF -- 11 xF-bins (-0.1-1.0), qT-differential, Q=[4.0,9.0] fixed
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/FNAL-537/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/piDY/"

### pi- beam on tungsten (W) target, modern nuclear code (per user instruction)
ps_def, h_1, h_2, proc_id = 2, 2, 184074, 1
### fixed 400 GeV/c pi- beam on a stationary nucleon target (as in old parsing)
s = 235.4
### 8% overall systematic (as in old parsing)
normErr = [0.08]
includeCuts = False


def _read_e537_table(path, tag_prefix):
    """
    Read an HEPData-style hand-typed table with comma-separated data rows,
    grouped into blocks by a "#: <tag>,,,<lo>TO<hi>" header line. Returns a
    list of (bin_low, bin_high, rows), rows = float lists [pT2_avg, pT2_low,
    pT2_high, xSec, err+, err-].
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
            label = line.split(",")[-1]
            bin_range = [float(v) for v in label.split("TO")]
            rows = []
        elif line == "" or line[0] in ("#", '"'):
            continue
        else:
            rows.append([float(v) for v in line.split(",")])
    if rows is not None:
        blocks.append((bin_range[0], bin_range[1], rows))
    return blocks


#%% -- Q-differential (E537(pi)-dQ)
blocks = _read_e537_table(path_to_data + "piminus+W(dQ).dat", "#: M(P=3 4) [GEV],,,")

ds = DataSet.empty("DY", name="E537.pi.dQ", comment="E537 pi-data Q-differential",
                    reference="Phys.Rev.D 93 (1988) 1377", normErr=normErr, isNormalized=False)

y_min, y_max = -0.1, 1.0  # full integrated xF acceptance, stashed in the y-fields
k = 0
for Q_min, Q_max, rows in blocks:
    for row in rows:
        qT_min, qT_max = row[1]**0.5, row[2]**0.5
        uncorrErr_0 = (row[4] - row[5]) / 2.
        if uncorrErr_0 <= 0:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"{ds.name}.{k}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            ### pb->nb unit conversion (as in old parsing)
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (Q_max - Q_min) * 0.001,
            includeCuts=includeCuts,
            xSec=row[3],
            uncorrErr_0=uncorrErr_0,
        ))
        k += 1

ds.save_csv(path_to_save + ds.name + ".csv")

#%% -- xF-differential (E537(pi)-dxF)
blocks = _read_e537_table(path_to_data + "pminus+W(dxF).dat", "#: XL(P=3 4,RF=CM,DEF=PL/PLMAX),,,")

ds = DataSet.empty("DY", name="E537.pi.dxF", comment="E537(pi) data xF-differential",
                    reference="Phys.Rev.D 93 (1988) 1377", normErr=normErr, isNormalized=False)

Q_min, Q_max = 4.0, 9.0  # fixed mass window
k = 0
for y_min, y_max, rows in blocks:
    for row in rows:
        qT_min, qT_max = row[1]**0.5, row[2]**0.5
        uncorrErr_0 = (row[4] - row[5]) / 2.
        if uncorrErr_0 <= 0:
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ds.add_point(dict(
            id=f"{ds.name}.{k}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            ### pb->nb unit conversion (as in old parsing)
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * 0.001,
            includeCuts=includeCuts,
            xSec=row[3],
            uncorrErr_0=uncorrErr_0,
        ))
        k += 1

ds.save_csv(path_to_save + ds.name + ".csv")
