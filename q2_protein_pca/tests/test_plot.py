# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json

import pandas as pd
from q2_protein_pca._plot import _generate_spec
from qiime2.plugin.testing import TestPluginBase

from q2_protein_pca.tests.data.expected_spec import (EXPECTED_SPEC,
                                                     EXPECTED_SPEC_WITH_NANS)


class PlotTests(TestPluginBase):

    package = 'q2_protein_pca.tests'

    def _prepare_data(self):
        data = [{'id': 'pos1', 'PC1': 0, 'PC2': 0,
                 'euclid_dist': 0, 'max_dist': 0.26, 'COL3': 1},
                {'id': 'pos2', 'PC1': 0.1, 'PC2': 0.2,
                 'euclid_dist': 0.05, 'max_dist': 0.26, 'COL3': 2},
                {'id': 'pos3', 'PC1': -0.5, 'PC2': -0.1,
                 'euclid_dist': 0.26, 'max_dist': 0.26, 'COL3': 3}]
        data_df = pd.DataFrame(data)
        sequence_ids = ["id1", "id2", "id3"]
        return data_df, sequence_ids

    def test_generate_spec(self):
        plot_data, sequence_ids = self._prepare_data()
        plot_data['COL3'] = plot_data['COL3'].astype(pd.Int64Dtype())

        obs_spec = _generate_spec(plot_data, 'PC1', 'PC2', sequence_ids)
        self.maxDiff = None
        self.assertDictEqual(obs_spec, EXPECTED_SPEC)
        json.dumps(obs_spec)

    def test_generate_spec_with_nans(self):
        plot_data, sequence_ids = self._prepare_data()
        plot_data.iloc[0, 5] = None

        obs_spec = _generate_spec(plot_data, 'PC1', 'PC2', sequence_ids)
        self.maxDiff = None
        self.assertDictEqual(obs_spec, EXPECTED_SPEC_WITH_NANS)
        json.dumps(obs_spec)
