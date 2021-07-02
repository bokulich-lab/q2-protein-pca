# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt
from q2_types.feature_data._transformer import AlignedProteinIterator

from q2_protein_pca._format import PositionMappingFormat
from q2_types.feature_data import AlignedProteinFASTAFormat, ProteinFASTAFormat
from qiime2.plugin.testing import TestPluginBase

from q2_protein_pca._alignment import map_positions, mafft


class AlignmentTests(TestPluginBase):

    package = 'q2_protein_pca.tests'

    def test_align_proteins(self):
        unaln_input_fp = self.get_data_path('raw-protein-sequences-3.fasta')
        unaln_input_seqs = ProteinFASTAFormat(unaln_input_fp, mode='r')

        aln_input_fp = self.get_data_path('aligned-protein-sequences-3.fasta')
        aln_input_seqs = AlignedProteinFASTAFormat(
            aln_input_fp, mode='r').view(AlignedProteinIterator)

        obs_seqs = mafft(unaln_input_seqs)

        for aln_seq, exp_seq in zip(
                obs_seqs.view(AlignedProteinIterator), aln_input_seqs):
            self.assertEqual(exp_seq.metadata['id'], aln_seq.metadata['id'])
            self.assertEqual(str(exp_seq), str(aln_seq))

    def test_map_positions(self):
        aln_input_fp = self.get_data_path('aligned-protein-sequences-2.fasta')
        aln_input_seqs = AlignedProteinFASTAFormat(
            aln_input_fp, mode='r').view(AlignedProteinIterator)

        exp_pos_fp = self.get_data_path('positions-mapping-2.csv')
        exp_pos = PositionMappingFormat(
            exp_pos_fp, mode='r').view(pd.DataFrame)

        obs_pos = map_positions(aln_input_seqs)

        pdt.assert_frame_equal(obs_pos, exp_pos)
