import pandas as pd
import pandas.testing as pdt
from q2_types.feature_data._transformer import AlignedProteinIterator

from q2_protein_pca._format import PositionMappingFormat
from q2_types.feature_data import AlignedProteinFASTAFormat
from qiime2.plugin.testing import TestPluginBase

from q2_protein_pca._alignment import map_positions


class AlignmentTests(TestPluginBase):

    package = 'q2_protein_pca.tests'

    def _prepare_sequences(self):
        aln_input_fp = self.get_data_path('aligned-protein-sequences-2.fasta')
        aln_input_seqs = AlignedProteinFASTAFormat(aln_input_fp, mode='r').view(AlignedProteinIterator)

        exp_pos_fp = self.get_data_path('positions-mapping-2.csv')
        exp_pos = PositionMappingFormat(exp_pos_fp, mode='r').view(pd.DataFrame)

        return aln_input_seqs, exp_pos

    def test_map_positions(self):
        aln_input_seqs, exp_pos = self._prepare_sequences()
        obs_pos = map_positions(aln_input_seqs)

        pdt.assert_frame_equal(obs_pos, exp_pos)
