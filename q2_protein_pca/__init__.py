# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._alignment import mafft, map_positions
from ._pca import pca
from ._plot import plot_loadings
from ._ranking import rank_alignment

__version__ = "2020.08"

__all__ = ['mafft', 'map_positions', 'pca', 'plot_loadings', 'rank_alignment']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
