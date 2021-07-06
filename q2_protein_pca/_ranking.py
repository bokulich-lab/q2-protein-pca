# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import aln_ranking as rank
import numpy as np
import pandas as pd
from q2_types.feature_data._transformer import AlignedProteinIterator

AA_MAP = {y: x for (x, y) in enumerate(list("-ABCDEFGHIKLMNPQRSTVWXYZ"))}


def _df_from_sequences(sequences: AlignedProteinIterator) -> pd.DataFrame:
    master_list, seq_ids = [], []
    for seq_record in sequences:
        master_list.append(list(seq_record.values.astype('str')))
        seq_ids.append(seq_record.metadata['id'])
    col_names = [f"pos{x+1}" for x in range(len(master_list[0]))]
    alignment_df = pd.DataFrame(master_list, columns=col_names, index=seq_ids)
    return alignment_df


def _get_occurrences(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(pd.value_counts).fillna(0).astype("int")


def _rank_columns(alignment_df: pd.DataFrame) -> pd.DataFrame:
    aln_unranked = alignment_df.replace(AA_MAP).astype(np.uint32).values
    aln_ranked = rank.rank_sequences(aln_unranked)
    aln_df_ranked = pd.DataFrame(
        aln_ranked, columns=alignment_df.columns, index=alignment_df.index)
    aln_df_ranked.index.name = "Sequence ID"
    return aln_df_ranked.astype("int")


def _rank(sequences: AlignedProteinIterator) -> pd.DataFrame:
    alignment_df = _df_from_sequences(sequences)
    ranking_table = _rank_columns(alignment_df)
    return ranking_table.astype("int")


def rank_alignment(sequences: AlignedProteinIterator) -> pd.DataFrame:
    return _rank(sequences)
