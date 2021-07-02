# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.util.testing as pdt

from qiime2.plugin.testing import TestPluginBase
from q2_protein_pca._format import (
    RankedProteinAlignmentFormat, PositionMappingFormat)


class TestTransformers(TestPluginBase):

    package = 'q2_protein_pca.tests'

    def setUp(self):
        super().setUp()
        self.protein_seqs = pd.DataFrame(
            [{'pos1': 2, 'pos2': 4, 'pos3': 2, 'pos4': 5},
             {'pos1': 2, 'pos2': 3, 'pos3': 1, 'pos4': 4},
             {'pos1': 2, 'pos2': 4, 'pos3': 2, 'pos4': 1}],
            index=pd.Index(['seq0', 'seq1', 'seq2'], name='Sequence ID')
        )
        self.position_map = pd.DataFrame(
            [{'seq0': 0, 'seq1': None, 'seq2': 0},
             {'seq0': None, 'seq1': None, 'seq2': 1},
             {'seq0': 1, 'seq1': 2, 'seq2': 3},
             {'seq0': 1, 'seq1': None, 'seq2': 3}],
            index=pd.Index([0, 1, 2, 3], name='Alignment position')
        ).astype('Int64')

    def test_ranked_aln_format_to_dataframe(self):
        _, obs = self.transform_format(
            RankedProteinAlignmentFormat, pd.DataFrame, 'protein-ranks.csv')
        self.assertIsInstance(obs, pd.DataFrame)
        pdt.assert_frame_equal(obs, self.protein_seqs)

    def test_dataframe_to_ranked_aln(self):
        transformer = self.get_transformer(
            pd.DataFrame, RankedProteinAlignmentFormat)

        obs = transformer(self.protein_seqs)
        self.assertIsInstance(obs, RankedProteinAlignmentFormat)

        obs_lines = list(open(str(obs)))
        self.assertEqual(obs_lines[0], 'Sequence ID,pos1,pos2,pos3,pos4\n')
        self.assertEqual(obs_lines[1], 'seq0,2,4,2,5\n')
        self.assertEqual(obs_lines[2], 'seq1,2,3,1,4\n')
        self.assertEqual(obs_lines[3], 'seq2,2,4,2,1\n')

    def test_position_map_format_to_dataframe(self):
        _, obs = self.transform_format(
            PositionMappingFormat, pd.DataFrame, 'positions-mapping-3.csv')
        self.assertIsInstance(obs, pd.DataFrame)
        pdt.assert_frame_equal(obs, self.position_map)

    def test_dataframe_to_position_map(self):
        transformer = self.get_transformer(
            pd.DataFrame, PositionMappingFormat)

        obs = transformer(self.position_map)
        self.assertIsInstance(obs, PositionMappingFormat)

        obs_lines = list(open(str(obs)))
        self.assertEqual(obs_lines[0], 'Alignment position,seq0,seq1,seq2\n')
        self.assertEqual(obs_lines[1], '0,0,,0\n')
        self.assertEqual(obs_lines[2], '1,,,1\n')
        self.assertEqual(obs_lines[3], '2,1,2,3\n')
        self.assertEqual(obs_lines[4], '3,1,,3\n')
