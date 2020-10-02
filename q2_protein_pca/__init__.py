# ----------------------------------------------------------------------------
# Copyright (c) 2017--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._alignment import mafft
from ._pca import pca
from ._ranking import rank_alignment

__version__ = "0.0.0-dev"

__all__ = ['mafft', 'pca', 'rank_alignment']
