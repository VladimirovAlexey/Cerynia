# Cerynia

*Cerynia (or Ceryneia) was a region of ancient Greece that was sacred to Artemide (Artemis),
and it was where creatures of Artemis, like the Ceryneian Hind or the Golden Hind, dwelt.*

Python library for loading, processing, and fitting experimental data together with
[artemide](https://github.com/VladimirovAlexey/artemide-public)/`harpy`, for TMD
(transverse momentum dependent) phenomenology.

Cerynia takes care of the "data side" of a fit: reading experimental datasets, applying
kinematic cuts, building covariance matrices, matching theory to data, and computing chi2 —
while `harpy` (the Python interface to `artemide`) provides the theory predictions.

## Requirements

- Python 3, `numpy`, `pandas`, `scipy`
- `harpy`/`artemide`, only needed for the parts of Cerynia that actually compute theory
  (`Cerynia.harpyInterface`, `Cerynia.saveTMDGrid`, `Cerynia.aTMDeReplicaSet.set`).
  Cerynia's data-handling classes (`Point`, `DataSet`, `DataMultiSet`) work without it.

## Repository layout

```
Cerynia/                    the Python package itself
    Point.py                 point schema/validation for DY, SIDIS, G2, D2 processes
    DataSet.py                a single experimental dataset (points + covariance + chi2)
    DataMultiSet.py            a collection of DataSets, combined chi2/theory-matching
    aTMDeReplicaSet.py         NP-parameter replica sets for artemide/Snowflake models
    harpyInterface.py          xsec()/chi2()/print_chi2_table() via harpy
    saveTMDGrid.py             write TMD/TMDFF grids in the TMDlib text-grid format

DataLib/                    the data library: experimental datasets as .csv files,
                             grouped by category (DY, DY_W, DY_angular, piDY, SIDIS,
                             Sivers, wgt, G2, D2); see DataLib/DataLib_Overview.xlsx
                             for the full list of included experiments

Replicas/                   NP-parameter replica sets (fit results) in JSON format,
                             loadable via Cerynia.aTMDeReplicaSet.from_json

FittingPrograms/             analysis/fitting scripts built on top of Cerynia
    EXAMPLE/                   minimal, heavily-commented usage examples — **start here**
    <fit-name>/                individual fitting projects

OtherPrograms/
    DataParsing/               scripts that generate the .csv files under DataLib/
    ParsingOldReplicaFiles/    .rep -> JSON converter for old artemide/Snowflake replicas
```

## Getting started

The best entry point is `FittingPrograms/EXAMPLE/`:

- `LoadData_xSec_chi2.py` — load artemide (`harpy`), load a replica set, load data from
  `DataLib/`, apply cuts, prepare it, compute cross-sections and chi2.
- `CreateDataFile.py` — build a `DataSet` from scratch point-by-point (no file I/O),
  explains which point fields are mandatory vs. optional, what `thFactor` and the
  normalization/`isNormalized` options mean, and how to save the result to a `.csv` file.

## Data format

Each dataset is a `.csv` file with set-level metadata (name, reference, `processType`,
normalization info) written as leading `# key: value` comment lines, followed by a normal
csv table of points — one row per data point, one column per kinematic variable/error
source. See `Cerynia/Point.py` for the schema of each `processType`, or call
`Cerynia.Point.schema("DY")` (etc.) to print it.

## License

MIT — see `LICENSE`.
