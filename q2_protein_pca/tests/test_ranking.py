# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import islice
import pandas as pd
import pandas.testing as pdt
from q2_types.feature_data import AlignedProteinFASTAFormat
from q2_types.feature_data._transformer import AlignedProteinIterator
from qiime2.plugin.testing import TestPluginBase

from q2_protein_pca import rank_alignment
from q2_protein_pca._format import RankedProteinAlignmentFormat
from q2_protein_pca._ranking import (
    _get_occurrences, _df_from_sequences, _rank_columns)


class RankingTests(TestPluginBase):

    package = 'q2_protein_pca.tests'

    def _prepare_sequences(self):
        input_fp = self.get_data_path('aligned-protein-sequences-1.fasta')
        input_sequences = AlignedProteinFASTAFormat(
            input_fp, mode='r').view(AlignedProteinIterator)

        ranks_fp = self.get_data_path('aligned-protein-ranks-1.csv')
        expected_ranks = RankedProteinAlignmentFormat(
            ranks_fp, mode='r').view(pd.DataFrame)

        return input_sequences, expected_ranks

    def test_df_from_sequences(self):
        input_seqs, _ = self._prepare_sequences()
        input_seqs = islice(input_seqs, 2)

        obs_df = _df_from_sequences(input_seqs)
        exp_df = pd.DataFrame(
            {"pos1": ["A", "A"], "pos2": ["A", "C"], "pos3": ["A", "C"],
             "pos4": ["C", "D"], "pos5": ["F", "G"], "pos6": ["C", "C"],
             "pos7": ["C", "D"], "pos8": ["D", "E"], "pos9": ["A", "A"]},
            index=["seq0", "seq1"])
        pdt.assert_frame_equal(obs_df, exp_df)

    def test_get_occurences(self):
        input_df = pd.DataFrame({"pos1": ["A", "A", "A", "A"],
                                 "pos2": ["-", "B", "B", "D"],
                                 "pos3": ["A", "C", "A", "-"],
                                 "pos4": ["-", "D", "B", "-"]},
                                index=["seq0", "seq1", "seq2", "seq3"])
        obs_occurences = _get_occurrences(input_df)
        exp_occurences = pd.DataFrame({"pos1": [0, 4, 0, 0, 0],
                                       "pos2": [1, 0, 2, 0, 1],
                                       "pos3": [1, 2, 0, 1, 0],
                                       "pos4": [2, 0, 1, 0, 1]},
                                      index=["-", "A", "B", "C", "D"])
        pdt.assert_frame_equal(obs_occurences, exp_occurences)

    def test_rank_columns(self):
        input_seqs = pd.DataFrame({"pos1": ["A", "A", "A", "A"],
                                   "pos2": ["-", "B", "B", "D"],
                                   "pos3": ["A", "C", "A", "-"],
                                   "pos4": ["-", "D", "B", "-"]},
                                  index=["seq0", "seq1", "seq2", "seq3"])

        obs_ranks = _rank_columns(input_seqs)
        exp_ranks = pd.DataFrame({"pos1": [1, 1, 1, 1],
                                  "pos2": [0, 2, 2, 1],
                                  "pos3": [2, 1, 2, 0],
                                  "pos4": [0, 1, 2, 0]},
                                 index=["seq0", "seq1", "seq2", "seq3"])
        exp_ranks.index.name = "Sequence ID"
        pdt.assert_frame_equal(obs_ranks, exp_ranks)

    def test_ranking(self):
        input_seqs, exp_ranks = self._prepare_sequences()
        obs_ranks = rank_alignment(input_seqs)
        pdt.assert_frame_equal(obs_ranks, exp_ranks)
