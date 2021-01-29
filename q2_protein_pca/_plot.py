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


def _generate_spec(plot_values: pd.DataFrame,
                   x_col_name: str,
                   y_col_name: str,
                   sequence_ids: list) -> dict:
    spec = {
        '$schema': 'https://vega.github.io/schema/vega/v4.2.json',
        'width': 300,
        'height': 300,
        'data': [
            {'name': 'values',
             'values': plot_values.replace(
                 {np.nan: None}).to_dict(orient='records')},
        ],
        'scales': [
            {'name': 'xScale',
             'domain': {'data': 'values', 'field': x_col_name},
             'range': 'width'},
            {'name': 'yScale',
             'domain': {'data': 'values', 'field': y_col_name},
             'range': 'height'}],
        'axes': [
            {'scale': 'xScale', 'orient': 'bottom', 'title': x_col_name},
            {'scale': 'yScale', 'orient': 'left', 'title': y_col_name}],
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
            },
            {
                "name": "sequenceID",
                "description": "Sequence ID",
                "bind": {
                    "input": "select",
                    "options": sequence_ids,
                    "element": "#sequence-id-selector"
                }
            },
            {
                "name": "hideMissingPositions",
                "description": "Hide missing positions",
                "bind": {
                    "input": "checkbox",
                    "element": "#hide-positions-selector"
                },
            }
        ],
        'marks': [
            {'type': 'symbol',
             'from': {'data': 'values'},
             'encode': {
                 'hover': {
                     'fill': {'value': '#d62728'},
                     'opacity': {'value': 0.8}},
                 'enter': {
                     'x': {'scale': 'xScale', 'field': x_col_name},
                     'y': {'scale': 'yScale', 'field': y_col_name}},
                 'update': {
                     'fill': [
                         {
                             'test': "datum.euclid_dist / datum.max_distance "
                                     "<= (1 - conservationLevel / 100)",
                             'value': '#3182bd'
                         },
                         {'value': 'black'}
                     ],
                     'opacity': [
                         {
                             'test': "datum[sequenceID] == null "
                                     "&& hideMissingPositions",
                             'value': 0.0
                         },
                         {'value': 0.8}],
                     'tooltip': {
                         'signal': f"{{'title': 'position ' + datum['id'], "
                                   f"'{x_col_name}': datum['{x_col_name}'], "
                                   f"'{y_col_name}': datum['{y_col_name}']}}"}
                 }}}]}
    return spec


def _plot_loadings(
        output_dir: str,
        pca_loadings_df: pd.DataFrame,
        positions_mapping: pd.DataFrame,
        pdb_id: str,
        nterm_offset: int):
    context = dict()

    # convert to 1-based indexing
    positions_mapping += 1
    context['position_data'] = positions_mapping.to_json(orient='records')

    positions_mapping.index = pca_loadings_df.index

    plot_values = pca_loadings_df.iloc[:, :2]
    plot_values.columns = [f"PC{x+1}" for x in plot_values.columns]
    x_col, y_col = 'PC1', 'PC2'

    plot_values.index.name = 'id'

    # calculate euclidean distances from (0,0) for each pair
    plot_values['euclid_dist'] = plot_values.apply(
        lambda x: _euclidean_distance(x[x_col], x[y_col]), axis=1)
    plot_values['max_distance'] = plot_values['euclid_dist'].max()

    plot_values = pd.concat([plot_values, positions_mapping], axis=1)
    plot_values.to_csv(os.path.join(output_dir, 'data.tsv'),
                       header=True, index=True, sep='\t')
    plot_values = plot_values.reset_index(drop=False)

    spec = _generate_spec(
        plot_values=plot_values,
        x_col_name=x_col,
        y_col_name=y_col,
        sequence_ids=list(
            positions_mapping.columns))

    context['vega_spec'] = json.dumps(spec)
    context['max_count'] = plot_values.shape[0]
    context['max_distance'] = plot_values['euclid_dist'].max()
    if pdb_id:
        context['pdb_id'] = pdb_id
        context['nterm_offset'] = nterm_offset

    copy_tree(os.path.join(TEMPLATES, 'loadings'), output_dir)

    index = os.path.join(TEMPLATES, 'loadings', 'index.html')
    q2templates.render(index, output_dir, context=context)


def plot_loadings(
        output_dir: str,
        pca_loadings: OrdinationResults,
        positions_mapping: pd.DataFrame,
        pdb_id: str = None,
        nterm_offset: int = 1) -> None:
    loadings_df = pca_loadings.samples
    _plot_loadings(
        output_dir, loadings_df, positions_mapping, pdb_id, nterm_offset)
