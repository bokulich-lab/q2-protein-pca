# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_types.feature_data._transformer import AlignedProteinIterator


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


def _rank_alphabet(aa_series: pd.Series) -> dict:
    ranks_df = pd.DataFrame({"aa": aa_series.index, "count": aa_series.values})

    # remove gaps and aas with no occurence at a given position
    ranks_df = ranks_df[(ranks_df["aa"] != "-") & (ranks_df["count"] != 0)]

    # sort by occurence and reverse-alphabetically
    ranks_df.sort_values(
        by=["count", "aa"], ascending=[True, False], inplace=True)

    # rank the aas according to the order
    ranks_df["rank"] = range(1, ranks_df.shape[0] + 1)

    # convert to dict
    ranks_dict = {row["aa"]: row["rank"] for index, row in ranks_df.iterrows()}

    # replace gap (-) counts with 0
    ranks_dict.update({"-": 0})

    return ranks_dict


def _rank_columns(alignment_df: pd.DataFrame) -> pd.DataFrame:
    occurence_freqs = _get_occurrences(alignment_df)

    # calculate AA ranks
    aa_ranks = occurence_freqs.apply(_rank_alphabet, axis=0)

    # replace original AAs with their ranks
    aa_ranks_dict = aa_ranks.to_dict()

    # aa_ranks_new = dict()
    # for key in aa_ranks_dict.keys():
    #     aa_ranks_new[f"pos{int(key)}"] = aa_ranks_dict[key]

    aln_df_ranked = alignment_df.replace(aa_ranks_dict)
    aln_df_ranked.index.name = "Sequence ID"
    return aln_df_ranked.astype("int")


def _rank(sequences: AlignedProteinIterator) -> pd.DataFrame:
    alignment_df = _df_from_sequences(sequences)
    ranking_table = _rank_columns(alignment_df)
    return ranking_table.astype("int")


def rank_alignment(sequences: AlignedProteinIterator) -> pd.DataFrame:
    return _rank(sequences)
