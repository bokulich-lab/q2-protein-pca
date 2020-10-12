import pandas as pd
import pandas.testing as pdt
from q2_types.feature_data import AlignedProteinFASTAFormat
from q2_types.feature_data._transformer import AlignedProteinIterator
from qiime2.plugin.testing import TestPluginBase

from q2_protein_pca import rank_alignment
from q2_protein_pca._format import RankedProteinAlignmentFormat


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

    def test_ranking(self):
        input_seqs, exp_ranks = self._prepare_sequences()
        obs_ranks = rank_alignment(input_seqs)

        pdt.assert_frame_equal(obs_ranks, exp_ranks)
