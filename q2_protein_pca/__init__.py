# ----------------------------------------------------------------------------
# Copyright (c) 2017--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._alignment import mafft, map_positions
from ._pca import pca
from ._plot import plot_loadings
from ._ranking import rank_alignment

__version__ = "2020.08"

__all__ = ['mafft', 'map_positions', 'pca', 'plot_loadings', 'rank_alignment']
