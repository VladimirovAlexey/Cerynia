"""
Cerynia — data processing library for TMD phenomenology fits.

Modules:
    Point           : schema definitions and validation for data points
    DataSet         : single-experiment dataset with chi2, systematic shifts, replica generation
    DataMultiSet    : collection of DataSets with combined chi2 and flat-vector interface
    aTMDeReplicaSet : NP-parameter replica sets for artemide/Snowflake (harpy) models
    saveTMDGrid     : write TMD/TMDFF grids to the standard TMDlib-style text format
"""

from . import Point
from .DataSet         import DataSet
from .DataMultiSet    import DataMultiSet
from .aTMDeReplicaSet import aTMDeReplicaSet

try:
    from . import harpyInterface
    from . import saveTMDGrid
except ImportError:
    pass  # harpy not installed; harpyInterface/saveTMDGrid unavailable
