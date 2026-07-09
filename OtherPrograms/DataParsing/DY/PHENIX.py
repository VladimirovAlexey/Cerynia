"""
Parsing of PHENIX 200 GeV (p+p) Drell-Yan invariant cross-section data.

Source: /data/arTeMiDe_Repository/data/Phenix/fig33.txt
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/PHE200.csv

Ported from the "PHENIX 200" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles(RHIC).py -- definitions
(process code, s, Q/y-range, thFactor, error split) are kept identical to the
old parsing.

NOTE: thFactor uses the same fixed-target-style "invariant cross section"
Jacobian as E288/E772/E605/E537 (1/(qT_max**2-qT_min**2)/(y_max-y_min)*
0.3183098861838*0.001, here with the *0.001 nb->pb factor matching the old
code's own comment, unlike E288's mismatched comment), even though this is a
p+p collider measurement (ps_def/h_1/h_2 = 1,1,1, not a fixed-target h-code).
Per the standing atmdeFactor rule, atmdeFactor is applied by multiplying the
old expression, same as every other DY case.

Raw format: whitespace-separated, 2 header lines (title + dashes), 12 data
rows, columns [pt, value, staterror, syserror(b,high), syserror(b,low),
syserror(c)].

Dataset:
    PHE200 -- arXiv:1805.02448
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"

#%%
# ============================================================================
# PHENIX 200 GeV -- arXiv:1805.02448
# ============================================================================
with open(path_to_data + "Phenix/fig33.txt") as f:
    lines = [line.rstrip("\n") for line in f]
rows = [[float(v) for v in line.split()] for line in lines[2:]]

### fixed low-mass window (off Z-peak, no Q_avg override, as in old parsing)
Q_min, Q_max = 4.8, 8.2
### Process is DY for p-p
ps_def, h_1, h_2, proc_id = 1, 1, 1, 1
### single forward-rapidity bin (as in old parsing)
y_min, y_max = 1.2, 2.2
### 12% normalization uncertainty (as in old parsing)
normErr = [0.12]

ds = DataSet.empty("DY", name="PHE200", comment="PHENIX 200GeV data", reference="arXiv:1805.02448",
                    normErr=normErr, isNormalized=False)

for i, (pt, xSec, staterr, sysb_hi, sysb_lo, sysc) in enumerate(rows):
    qT_min, qT_max = pt - 0.25, pt + 0.25  ### bin width injected manually (table gives bin center only)

    #### artemide computes the avarage over bin, this factor compensates this
    atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

    ds.add_point(dict(
        id=f"PHE200.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
        s=200.**2, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
        qT_min=qT_min, qT_max=qT_max,
        ### invariant cross-section Jacobian (as in old parsing); 0.3183098861838 = 1/pi, 0.001 = nb->pb
        thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * 0.3183098861838 * 0.001,
        includeCuts=False,
        xSec=xSec,
        uncorrErr_0=staterr,
        ### "b"-flavor systematic, both magnitudes positive; average them
        uncorrErr_1=(sysb_hi + sysb_lo) / 2.,
        ### "c"-flavor systematic, treated as correlated (as in old parsing)
        corrErr_0=sysc,
    ))

ds.save_csv(path_to_save + ds.name + ".csv")
