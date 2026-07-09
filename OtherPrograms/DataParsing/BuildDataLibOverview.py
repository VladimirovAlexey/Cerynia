"""
Builds a human-readable Excel overview of every parsed dataset in Cerynia/DataLib,
grouped by DataLib subfolder/category (e.g. DY, DY_W, ...), then into
subsections by experiment/source file within each category.

Output: /data/arTeMiDe_Repository/Cerynia/DataLib/DataLib_Overview.xlsx

Re-run this whenever new datasets are added to DataLib/ -- it re-derives all
numeric columns (point count, Q-range, pT-range, process codes) directly from
the saved CSVs, so it never goes stale except for the hand-written columns
(experiment name, physical-process string, comments), which are looked up
from the <CATEGORY>_SECTIONS tables below and must be extended by hand for
new datasets. To add a new category (new DataLib subfolder), add a
"<CATEGORY>_SECTIONS" list following the same pattern as DY_SECTIONS /
DY_W_SECTIONS, then add it to the CATEGORIES list with its subfolder name.
"""

import sys
sys.path.append("/data/arTeMiDe_Repository/Cerynia")

from Cerynia import DataSet
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

DATALIB_ROOT = "/data/arTeMiDe_Repository/Cerynia/DataLib/"
OUTPUT = "/data/arTeMiDe_Repository/Cerynia/DataLib/DataLib_Overview.xlsx"

COLUMNS = [
    "Dataset", "Experiment", "Reference", "Process", "Process (artemide)",
    "Normalized", "N points", "Q min [GeV]", "Q max [GeV]",
    "pT/qT min [GeV]", "pT/qT max [GeV]", "Comments",
]

# Each category: (category_title, subfolder_under_DataLib, [section, ...])
# Each section: (title, [(dataset_name, experiment, physical_process, comment), ...])
# comment is "" unless there is a known/open issue to flag.
DY_SECTIONS = [
    ("Tevatron (CDF / D0)", [
        ("CDF1", "CDF",  "p + pbar -> Z/gamma* -> l+l-", ""),
        ("CDF2", "CDF",  "p + pbar -> Z/gamma* -> l+l-", ""),
        ("D01",  "D0",   "p + pbar -> Z/gamma* -> l+l-", ""),
        ("D02",  "D0",   "p + pbar -> Z/gamma* -> l+l-", ""),
        ("D02m", "D0",   "p + pbar -> Z/gamma* -> l+l-", ""),
    ]),
    ("ATLAS 7 TeV", [
        ("A7-00y10", "ATLAS (7 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A7-10y20", "ATLAS (7 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A7-20y24", "ATLAS (7 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
    ]),
    ("ATLAS 8 TeV", [
        ("A8-00y04",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-04y08",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-08y12",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-12y16",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-16y20",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-20y24",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-46Q66",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-116Q150", "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
    ]),
    ("ATLAS 8 TeV, normalized to 1/sigma", [
        ("A8-00y04-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-04y08-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-08y12-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-12y16-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-16y20-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-20y24-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-46Q66-norm",   "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("A8-116Q150-norm", "ATLAS (8 TeV)", "p + p -> Z/gamma* -> l+l-",
         "same underlying measurement as A8-* (absolute), normalized to 1/sigma instead"),
    ]),
    ("ATLAS 13 TeV, normalized to 1/sigma", [
        ("A13-norm", "ATLAS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
    ]),
    ("CMS", [
        ("CMS7", "CMS (7 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS8", "CMS (8 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
    ]),
    ("CMS 13 TeV, y-differential", [
        ("CMS13-00y04", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-04y08", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-08y12", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-12y16", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-16y24", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-00y04-norm", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-04y08-norm", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-08y12-norm", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-12y16-norm", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13-16y24-norm", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-",
         "pT bin edges restored from a hardcoded table (raw file only gives bin-center labels), not re-derivable from source"),
    ]),
    ("CMS 13 TeV, large Q (preliminary)", [
        ("CMS13_dQ_50to76",    "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13_dQ_76to106",   "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13_dQ_106to170",  "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13_dQ_170to350",  "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-", ""),
        ("CMS13_dQ_350to1000", "CMS (13 TeV)", "p + p -> Z/gamma* -> l+l-",
         "preliminary CMS-PAS result (hand-typed table); thFactor uses hardcoded per-point FSR correction factors, not re-derivable from source"),
    ]),
    ("LHCb", [
        ("LHCb7",  "LHCb (7 TeV)",  "p + p -> Z/gamma* -> l+l-", ""),
        ("LHCb8",  "LHCb (8 TeV)",  "p + p -> Z/gamma* -> l+l-", ""),
        ("LHCb13(2016)", "LHCb (13 TeV)", "p + p -> Z/gamma* -> l+l-",
         "hand-typed table (pre-official HepData release); provisional. Renamed from plain 'LHCb13' to disambiguate from LHCb13(2021)"),
    ]),
    ("LHCb 13 TeV, 2021 update", [
        ("LHCb13_dy", "LHCb (13 TeV)", "p + p -> Z/gamma* -> l+l-",
         "double-differential in (y,qT); '(2021)' label dropped from the name per user instruction"),
        ("LHCb13(2021)", "LHCb (13 TeV)", "p + p -> Z/gamma* -> l+l-",
         "y-integrated companion to LHCb13_dy (same paper); '(2021)' label kept per user instruction"),
    ]),
    ("RHIC (PHENIX / STAR)", [
        ("PHE200", "PHENIX (200 GeV)", "p + p -> gamma* -> l+l-",
         "invariant-xsec thFactor (qT^2 Jacobian, fixed-target-style formula) despite p+p collider process code"),
        ("STAR510", "STAR (510 GeV)", "p + p -> Z/gamma* -> l+l-",
         "qT bin edges restored from a hardcoded table (raw file only gives measured centroid), not re-derivable from source"),
    ]),
    ("E288 (fixed-target)", [
        ("E288-200", "E288", "p + Cu -> gamma* -> mu+mu- X",
         "ps_def corrected 2->1 per DataProcessor.git 'Corrections in E228' (2026-01-02). thFactor comment still says '0.001 for nb' but *1000 is applied; ported as coded, not corrected"),
        ("E288-300", "E288", "p + Cu -> gamma* -> mu+mu- X",
         "ps_def corrected 2->1 (same fix as E288-200). thFactor comment mismatch (see E288-200) still unresolved"),
        ("E288-400", "E288", "p + Cu -> gamma* -> mu+mu- X",
         "ps_def corrected 2->1 (same fix as E288-200). thFactor comment mismatch (see E288-200) still unresolved"),
    ]),
    ("E772 (fixed-target)", [
        ("E772", "E772", "p + D -> gamma* -> mu+mu- X",
         "thFactor gained the missing Feynman-x Jacobian (ffactor) per DataProcessor.git 'Corrected thFactor in E605 and E772' (2026-01-05)"),
    ]),
    ("E605 (fixed-target)", [
        ("E605", "E605", "p + Cu -> gamma* -> mu+mu- X",
         "old file mislabeled this target with proc_id=102 (the E772/deuteron code); fixed. thFactor also gained the missing Feynman-x Jacobian (ffactor, xF=0.1 fixed), same 2026-01-05 commit as E772"),
    ]),
    ("E537 (fixed-target)", [
        ("E537-dQ",  "E537", "pbar + W -> gamma* -> mu+mu- X", ""),
        ("E537-dxF", "E537", "pbar + W -> gamma* -> mu+mu- X",
         "h_1/h_2 order was swapped vs dQ in the original main-file parsing; fixed to match dQ (h_1=-1 pbar, h_2=184074 W), per DataProcessor.git's low-energy.py using that order consistently"),
    ]),
]

DY_W_SECTIONS = [
    ("Tevatron (D0 / CDF), run 1", [
        ("D01_W", "D0 (1.8 TeV)", "p + pbar -> W -> e nu",
         "Q-window is qT-dependent (kinematic bound M^2>E1+E2, per old parsing's comment), not a fixed range. Renamed from old 'D0run1-W' to match DY branch's D01 naming (same experiment, W final state)"),
        ("CDF1_W", "CDF (1.8 TeV)", "p + pbar -> W -> e nu",
         "qT bins restored from a hardcoded variable-width binning function around each point's centroid ('high-qT bins assumed from the plot'), not re-derivable from source. Rapidity restriction placed in cutParams (not y, unlike D01_W) -- inconsistent between the two datasets but ported as-is. Renamed from old 'CDFrun1-W' to match DY branch's CDF1 naming"),
    ]),
    ("CMS 8 TeV", [
        ("CMS8_W-electron", "CMS (8 TeV)", "p + p -> W -> e nu", ""),
        ("CMS8_W-muon", "CMS (8 TeV)", "p + p -> W -> mu nu",
         "cutParams' eta window (-2.5,2.5) does not match this dataset's own y-window (-2.1,2.1); ported verbatim from old parsing, not corrected"),
    ]),
]

SIDIS_SECTIONS = [
    ("HERMES 3D (z,x,pT), proton target", [
        ("hermes3D.p.pi+", "HERMES", "e + p -> e + pi+ X", ""),
        ("hermes3D.p.pi-", "HERMES", "e + p -> e + pi- X", ""),
        ("hermes3D.p.k+",  "HERMES", "e + p -> e + K+ X", ""),
        ("hermes3D.p.k-",  "HERMES", "e + p -> e + K- X", ""),
        ("hermes3D.p.pi+.no-vmsub", "HERMES", "e + p -> e + pi+ X", "no vector-meson-background subtraction (vmsub is the default/unlabeled case above)"),
        ("hermes3D.p.pi-.no-vmsub", "HERMES", "e + p -> e + pi- X", "no vector-meson-background subtraction"),
        ("hermes3D.p.k+.no-vmsub",  "HERMES", "e + p -> e + K+ X", "no vector-meson-background subtraction"),
        ("hermes3D.p.k-.no-vmsub",  "HERMES", "e + p -> e + K- X",
         "no vector-meson-background subtraction. All hermes3D sets: rows with all-zero fields (unpopulated (x,z,pT) grid cells) dropped per user confirmation; thFactor read from a raw DIS-normalization column (MSHT20nnlo), not computed from bin widths alone"),
    ]),
    ("HERMES 3D (z,x,pT), deuteron target", [
        ("hermes3D.d.pi+", "HERMES", "e + d -> e + pi+ X", ""),
        ("hermes3D.d.pi-", "HERMES", "e + d -> e + pi- X", ""),
        ("hermes3D.d.k+",  "HERMES", "e + d -> e + K+ X", ""),
        ("hermes3D.d.k-",  "HERMES", "e + d -> e + K- X", ""),
        ("hermes3D.d.pi+.no-vmsub", "HERMES", "e + d -> e + pi+ X", "no vector-meson-background subtraction"),
        ("hermes3D.d.pi-.no-vmsub", "HERMES", "e + d -> e + pi- X", "no vector-meson-background subtraction"),
        ("hermes3D.d.k+.no-vmsub",  "HERMES", "e + d -> e + K+ X", "no vector-meson-background subtraction"),
        ("hermes3D.d.k-.no-vmsub",  "HERMES", "e + d -> e + K- X", "no vector-meson-background subtraction. M_target uses the proton mass even for this deuteron-target set (per-nucleon SIDIS convention, ported verbatim)"),
    ]),
    ("COMPASS, deuteron target", [
        ("compass.d.h+", "COMPASS", "mu + d -> mu + h+ X", ""),
        ("compass.d.h-", "COMPASS", "mu + d -> mu + h- X",
         "isoscalar unidentified charged hadron; M_product uses the pion mass as a stand-in (as in old parsing)"),
    ]),
]

CATEGORIES = [
    ("DY -- unpolarized Drell-Yan", "DY", DY_SECTIONS),
    ("DY_W -- W-boson production", "DY_W", DY_W_SECTIONS),
    ("SIDIS -- unpolarized SIDIS", "SIDIS", SIDIS_SECTIONS),
]

# -- Build workbook ----------------------------------------------------------

wb = Workbook()
ws = wb.active
ws.title = "DataLib overview"

title_font     = Font(bold=True, size=14)
category_font  = Font(bold=True, size=13, color="FFFFFF")
category_fill  = PatternFill("solid", fgColor="203864")
section_font   = Font(bold=True, size=11, color="FFFFFF")
section_fill   = PatternFill("solid", fgColor="4472C4")
header_font    = Font(bold=True)
header_fill    = PatternFill("solid", fgColor="D9E1F2")
thin           = Side(style="thin", color="BFBFBF")
border         = Border(left=thin, right=thin, top=thin, bottom=thin)
wrap           = Alignment(wrap_text=True, vertical="top")
center         = Alignment(horizontal="center", vertical="center")

ncols = len(COLUMNS)

ws.cell(row=1, column=1, value="Cerynia DataLib -- parsed datasets overview").font = title_font
row = 3

for category_title, subfolder, sections in CATEGORIES:
    # Category banner
    ws.cell(row=row, column=1, value=category_title)
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=ncols)
    cell = ws.cell(row=row, column=1)
    cell.font = category_font
    cell.fill = category_fill
    cell.alignment = Alignment(horizontal="left", vertical="center")
    row += 2

    datalib_dir = DATALIB_ROOT + subfolder + "/"

    for section_title, entries in sections:
        # Section banner
        ws.cell(row=row, column=1, value=section_title)
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=ncols)
        cell = ws.cell(row=row, column=1)
        cell.font = section_font
        cell.fill = section_fill
        cell.alignment = Alignment(horizontal="left", vertical="center")
        row += 1

        # Header row
        for c, name in enumerate(COLUMNS, start=1):
            cell = ws.cell(row=row, column=c, value=name)
            cell.font = header_font
            cell.fill = header_fill
            cell.border = border
            cell.alignment = center
        row += 1

        for dsname, experiment, process, comment in entries:
            ds = DataSet.from_csv(datalib_dir + dsname + ".csv")
            df = ds.df
            procs = df[["ps_def", "h_1", "h_2", "proc_id"]].drop_duplicates().values.tolist()
            proc_str = "; ".join(str(p) for p in procs)
            ### DY branch uses qT_min/qT_max; SIDIS branch uses pT_min/pT_max instead
            pt_min_col = "qT_min" if "qT_min" in df.columns else "pT_min"
            pt_max_col = "qT_max" if "qT_max" in df.columns else "pT_max"

            values = [
                dsname,
                experiment,
                ds.reference,
                process,
                proc_str,
                "Yes" if ds.isNormalized else "",
                len(df),
                round(float(df["Q_min"].min()), 4),
                round(float(df["Q_max"].max()), 4),
                round(float(df[pt_min_col].min()), 4),
                round(float(df[pt_max_col].max()), 4),
                comment,
            ]
            for c, v in enumerate(values, start=1):
                cell = ws.cell(row=row, column=c, value=v)
                cell.border = border
                if c in (7, 8, 9, 10, 11):
                    cell.alignment = center
                elif c == 12:
                    cell.alignment = wrap
            row += 1

        row += 1  # blank separator row between sections

    row += 1  # extra blank row between categories

# -- Column widths, freeze header ------------------------------------------

widths = [14, 16, 26, 30, 20, 11, 9, 11, 11, 11, 11, 55]
for c, w in enumerate(widths, start=1):
    ws.column_dimensions[get_column_letter(c)].width = w

ws.freeze_panes = "A4"

wb.save(OUTPUT)
print(f"Wrote {OUTPUT}")
