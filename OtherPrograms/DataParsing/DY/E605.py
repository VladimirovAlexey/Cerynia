"""
Parsing of E605 (fixed-target, mu+mu- Drell-Yan) invariant cross-section data.

Source: /data/arTeMiDe_Repository/data/E605/E605_78.dat, E605_89.dat,
        E605_1011.dat, E605_1113.dat, E605_1318.dat  (YODA Scatter2D tables)
Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DY/E605.csv

Ported from the "E605" block of the old
DataProcessor/OtherPrograms/DataParsing/ReadDYdataFiles.py -- definitions
(process code, s, Q-range, y-range, thFactor, error split) are kept identical
to the old parsing. All 5 raw files feed into the single "E605" dataset, one
fixed Q-window per file, as in old parsing.

NOTE: thFactor is the fixed-target "invariant cross section" Jacobian,
1/(qT_max**2-qT_min**2)/(y_max-y_min)*ffactor -- structurally different from
the collider bin-integrated k/(qT_max-qT_min) form. Per user confirmation,
atmdeFactor still applies here: the old expression is multiplied by
atmdeFactor=(Q_max-Q_min)*(y_max-y_min)*(qT_max-qT_min), same as every other
DY case.

CORRECTED (per audit of DataProcessor.git history): the originally-ported
ReadDYdataFiles.py block used a plain constant 0.3183098861838 (=1/pi) in
place of ffactor -- missing a Feynman-x Jacobian. Commit "Corrected thFactor
in E605 and E772" (2026-01-05, DataProcessor's own repo) fixed this in a
parallel script (ReadDYdataFiles(low-energy).py) that never got merged back
into ReadDYdataFiles.py -- same commit that fixed the proc_id=102 mistake
noted above. ffactor = sqrt(4*(Q_avg**2+qT_avg**2)/s + xF**2) *
0.3183098861838, where xF=0.1 is a fixed value "specified in the text" (i.e.
E605's own paper), not derived from y_min/y_max -- ported verbatim.

Raw format: YODA Scatter2D ("xval xerr- xerr+ yval yerr- yerr+", no leading
tab, "BEGIN/IsRef/Path/Title/Type/# xval" header (6 lines), "END YODA_SCATTER2D"
+ trailing blank footer (2 lines) -- same format as CMS_8.dat. xerr-/xerr+ are
always 0 in these files (no real bin-width info), so qT bounds are manually
injected as center +- 0.1 GeV, as old parsing did. yerr-/yerr+ are both
positive magnitudes (not signed +/-); old parsing summed and halved them.

Dataset:
    E605 -- Phys.Rev.D 43 (1991) 2815
"""

import sys
from math import sqrt
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet

path_to_data = "/data/arTeMiDe_Repository/data/"
path_to_save = "/data/arTeMiDe_Repository/Cerynia/DataLib/DY/"


def _read_yoda_table(path, n_header=6):
    """
    Read a YODA Scatter2D table: "xval xerr- xerr+ yval yerr- yerr+", tab
    separated, no leading tab. Drops the n_header header lines and the
    trailing "END YODA_SCATTER2D" + blank line footer.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    lines = lines[n_header:-2]
    return [[float(v) for v in line.split("\t")] for line in lines]


#%%
# ============================================================================
# E605 -- Phys.Rev.D 43 (1991) 2815
# ============================================================================
### 800 GeV fixed-target beam momentum -> sqrt(s)
s = 38.76**2
### Process is virtual-photon DY for p-copper (target treated as isoscalar proton)
####ps_def, h_1, h_2, proc_id = 2, 1, 1, 102  #### there was a mistake in the old file (102=p+d)  ##(<6.04 definition)
ps_def, h_1, h_2, proc_id = 2, 1, 63029, 1    ##(>6.04 definition)

y_min, y_max = -0.1, 0.2
### 15% normalization uncertainty (as in old parsing)
normErr = [0.15]
### flat 5% systematic, applied per point as sysError*xSec (as in old parsing)
sysError = 0.05
### Feynman-x, fixed value "specified in the text" (as in old parsing's fix, not y_min/y_max-derived)
xF = 0.1

ds = DataSet.empty("DY", name="E605", comment="E605 data", reference="Phys.Rev.D 43 (1991) 2815",
                    normErr=normErr, isNormalized=False)

for fname, Q_min, Q_max, tag in [
    ("E605_78.dat",   7.0,  8.0,  "7Q8"),
    ("E605_89.dat",   8.0,  9.0,  "8Q9"),
    ("E605_1011.dat", 10.5, 11.5, "10Q11"),
    ("E605_1113.dat", 11.5, 13.5, "11Q13"),
    ("E605_1318.dat", 13.5, 18.0, "13Q18"),
]:
    rows = _read_yoda_table(path_to_data + "E605/" + fname)
    for i, (xval, xerrm, xerrp, yval, yerrm, yerrp) in enumerate(rows):
        qT_min, qT_max = xval - 0.1, xval + 0.1  ### xerr-/xerr+ are always 0; bin width injected manually
        xSec = yval
        if xSec == -50:  # missing-data sentinel (as in old parsing; never triggers on this raw data)
            continue

        #### artemide computes the avarage over bin, this factor compensates this
        atmdeFactor = (Q_max - Q_min) * (y_max - y_min) * (qT_max - qT_min)

        ### invariant cross-section Jacobian, corrected for Feynman-x (as in the
        ### 2026-01-05 fix); 0.3183098861838 = 1/pi
        Q_avg = (Q_min + Q_max) * 0.5
        qT_avg = (qT_min + qT_max) * 0.5
        ffactor = sqrt(4. * (Q_avg**2 + qT_avg**2) / s + xF**2) * 0.3183098861838

        ds.add_point(dict(
            id=f"E605.{tag}.{i}", ps_def=ps_def, h_1=h_1, h_2=h_2, proc_id=proc_id,
            s=s, Q_min=Q_min, Q_max=Q_max, y_min=y_min, y_max=y_max,
            qT_min=qT_min, qT_max=qT_max,
            thFactor=atmdeFactor / (qT_max**2 - qT_min**2) / (y_max - y_min) * ffactor,
            includeCuts=False,
            xSec=xSec,
            ### yerr-/yerr+ are both positive magnitudes in this format; average them
            uncorrErr_0=(yerrm + yerrp) / 2.,
            uncorrErr_1=sysError * xSec,
        ))

ds.save_csv(path_to_save + ds.name + ".csv")
