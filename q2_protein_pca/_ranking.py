import pandas as pd
import skbio
from q2_types.feature_data import AlignedProteinFASTAFormat, RankedProteinAlignmentFormat


def _get_occurrences(df):
    return df.apply(pd.value_counts).fillna(0).astype("int64")


def _rank_alphabet(pd_series):
    # remove gaps beofore ranking
    pd_series_1 = pd_series.loc[[x for x in pd_series.index if x != "-"]]
    # sort AAs according to occurrence (descending)
    temp_df = pd.DataFrame(
        {"idx_col": pd_series_1.index, pd_series_1.name: pd_series_1}
    )
    test_rank_x = temp_df.sort_values(
        by=[pd_series_1.name, "idx_col"], ascending=[False, False]
    )
    test_rank = pd.Series(
        data=test_rank_x[pd_series_1.name].values,
        index=test_rank_x["idx_col"].values,
    )

    # group by count of occurrences
    all_occurrences = test_rank.groupby(test_rank).agg("count")
    # list how many AAs with given count (no of occur, no of aas) in a count descending order
    list_of_rep_aas = sorted(
        list(all_occurrences.iteritems()), key=lambda x: x[0], reverse=True
    )
    # generates a list of tuples containing the amino acid and its occurrence frequency, e.g.:
    # [('A', 10), ('C', 10), ('S', 0), ('R', 0), ('B', 0), ('D', 0) ... ]
    all_items = list(test_rank.iteritems())
    all_items_ranked = []
    all_ranks_by_sum = 0
    for i, e in enumerate(list_of_rep_aas):
        all_ranks_by_sum += e[1] if e[0] != 0 else False

    cnt1 = 0
    tot_aas = all_ranks_by_sum
    for x, rep_aa in enumerate(list(list_of_rep_aas)):
        cnt2 = 0
        rep_aa_count, rep_aa_aanum = rep_aa
        if rep_aa_count != 0.0:
            list_of_ranks = []
            for i in range(rep_aa_aanum):
                rank = tot_aas - rep_aa_aanum + cnt2 + 1
                list_of_ranks.append(rank)
                cnt2 += 1
            tot_aas = tot_aas - rep_aa_aanum

            for j, j_rank in enumerate(list_of_ranks):
                aas_with_rank = (all_items[cnt1][0], j_rank)
                all_items_ranked.append(aas_with_rank)
                cnt1 += 1
    ranks_dict = dict(all_items_ranked)

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
    return aln_df_ranked


def _rank(sequences: str) -> RankedProteinAlignmentFormat:
    # for seq in skbio.io.read(sequences, format='fasta',
    #                           constructor=skbio.Protein):
    #     id = seq.metadata['id']

    master_list, seq_ids = [], []
    for seq_record in skbio.io.registry.read(sequences, format='fasta'):
    # for seq_record in SeqIO.parse(sequences, 'fasta'):
        master_list.append(list(seq_record.values.astype('str')))
        seq_ids.append(seq_record.metadata['id'])
    col_names = [f"pos{x+1}" for x in range(len(master_list[0]))]
    alignment_df = pd.DataFrame(master_list, columns=col_names, index=seq_ids)
    ranking_table = _rank_columns(alignment_df)

    result = RankedProteinAlignmentFormat()
    result_fp = str(result)

    ranking_table.to_csv(result_fp, sep=",", header=True, index=True)

    return result


def rank_alignment(sequences: AlignedProteinFASTAFormat) -> RankedProteinAlignmentFormat:
    sequences_fp = str(sequences)
    return _rank(sequences_fp)