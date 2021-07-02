# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_protein_pca._plot import _generate_spec
from qiime2.plugin.testing import TestPluginBase

from q2_protein_pca.tests.data.expected_spec import EXPECTED_SPEC


class PlotTests(TestPluginBase):

    package = 'q2_protein_pca.tests'

    def _prepare_data(self):
        data = [{'PC1': 0, 'PC2': 0},
                {'PC1': 0.1, 'PC2': 0.2},
                {'PC1': -0.5, 'PC2': -0.1}]
        data_df = pd.DataFrame(data, index=["pos1", "pos2", "pos3"])
        sequence_ids = ["id1", "id2", "id3"]
        return data_df, sequence_ids

    def test_generate_spec(self):
        plot_data, sequence_ids = self._prepare_data()

        obs_spec = _generate_spec(plot_data, 'PC1', 'PC2', sequence_ids)
        self.maxDiff = None
        self.assertDictEqual(obs_spec, EXPECTED_SPEC)
