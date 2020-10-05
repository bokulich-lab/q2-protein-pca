import json
import numpy as np
import os
from distutils.dir_util import copy_tree

import pkg_resources
import q2templates
from skbio import OrdinationResults
import pandas as pd

TEMPLATES = pkg_resources.resource_filename('q2_protein_pca', 'assets')


def _euclidean_distance(x, y):
    return np.sqrt(x**2 + y**2)

def _generate_spec(output_dir: str, pca_loadings_df: pd.DataFrame):
    plot_values = pca_loadings_df.iloc[:, :2]
    plot_values.columns = [f"PC{x+1}" for x in plot_values.columns]
    x_col, y_col = 'PC1', 'PC2'

    plot_values.index.name = 'id'

    # calculate euclidean distances from (0,0) for each pair
    plot_values['euclid_dist'] = plot_values.apply(lambda x: _euclidean_distance(x[x_col], x[y_col]), axis=1)
    plot_values['max_distance'] = plot_values['euclid_dist'].max()
    plot_values.to_csv(os.path.join(output_dir, 'data.tsv'),
                           header=True, index=True, sep='\t')
    plot_values = plot_values.reset_index(drop=False)

    spec = {
        '$schema': 'https://vega.github.io/schema/vega/v5.json',
        'width': 250,
        'height': 250,
        'data': [
            {'name': 'values',
             'values': plot_values.to_dict(orient='records')},
        ],
        'scales': [
            {'name': 'xScale',
             'domain': {'data': 'values', 'field': x_col},
             'range': 'width'},
            {'name': 'yScale',
             'domain': {'data': 'values', 'field': y_col},
             'range': 'height'}],
        'axes': [
            {'scale': 'xScale', 'orient': 'bottom', 'title': x_col},
            {'scale': 'yScale', 'orient': 'left', 'title': y_col}],
        'signals': [
            {
            "name": "conservationLevel",
            "description": "Conservation level [%]",
            "value": 90,
            "bind": {
                "input": "range",
                "min": 0,
                "max": 100,
                "step": 1,
                "debounce": 10,
                "element": "#conservation-level-slider"
                }
            }
        ],
        'marks': [
            {'type': 'symbol',
             'from': {'data': 'values'},
             'encode': {
                 'hover': {
                     'fill': {'value': '#BD0909'},
                     'opacity': {'value': 0.55}},
                 'enter': {
                     'x': {'scale': 'xScale', 'field': x_col},
                     'y': {'scale': 'yScale', 'field': y_col}},
                 'update': {
                     'fill': [
                         {
                            'test': "datum.euclid_dist / datum.max_distance <= (1 - conservationLevel / 100)",
                            'value': 'red'
                         },
                         {'value': 'black'}
                    ],
                     'opacity': {'value': 0.55},
                     'tooltip': {
                         'signal': "{{'title': 'position ' + datum['id'], '{0}': "
                                   "datum['{0}'], '{1}': datum['{1}']}}".format(x_col, y_col)}
                 }}}]}
    context = dict()
    context['vega_spec'] = json.dumps(spec)
    context['max_count'] = plot_values.shape[0]
    context['max_distance'] = plot_values['euclid_dist'].max()
    context['pca_data'] = plot_values.to_json(orient='records')

    copy_tree(os.path.join(TEMPLATES, 'loadings'), output_dir)

    index = os.path.join(TEMPLATES, 'loadings', 'index.html')
    q2templates.render(index, output_dir, context=context)


def plot_loadings(output_dir: str, pca_loadings: OrdinationResults) -> None:
    loadings_df = pca_loadings.samples
    _generate_spec(output_dir, loadings_df)
