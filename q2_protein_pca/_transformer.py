# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from ._format import RankedProteinAlignmentFormat, PositionMappingFormat
from q2_protein_pca.plugin_setup import plugin


@plugin.register_transformer
def _1(ff: RankedProteinAlignmentFormat) -> pd.DataFrame:
    with ff.open() as fh:
        df = pd.read_csv(fh, sep=",", index_col="Sequence ID")
        df = df.astype('int64')
        return df


@plugin.register_transformer
def _2(data: pd.DataFrame) -> RankedProteinAlignmentFormat:
    ff = RankedProteinAlignmentFormat()
    with ff.open() as fh:
        data = data.astype('int64')
        data.to_csv(fh, sep=",", header=True, index=True)
        return ff


@plugin.register_transformer
def _3(ff: PositionMappingFormat) -> pd.DataFrame:
    with ff.open() as fh:
        df = pd.read_csv(fh, sep=",", index_col='Alignment position')
        df = df.astype('Int64')
        return df


@plugin.register_transformer
def _4(data: pd.DataFrame) -> PositionMappingFormat:
    ff = PositionMappingFormat()
    with ff.open() as fh:
        data = data.astype('Int64')
        data.to_csv(fh, sep=",", header=True, index=True)
        return ff
