"""
Cerynia — data processing library for TMD phenomenology fits.

Modules:
    Point       : schema definitions and validation for data points
    DataSet     : single-experiment dataset with chi2, systematic shifts, replica generation
    DataMultiSet: collection of DataSets with combined chi2 and flat-vector interface
"""

from . import Point
from .DataSet      import DataSet
from .DataMultiSet import DataMultiSet

try:
    from . import harpyInterface
except ImportError:
    pass  # harpy not installed; harpyInterface unavailable
