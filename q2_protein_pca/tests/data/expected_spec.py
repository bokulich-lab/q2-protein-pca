# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

EXPECTED_SPEC = {
    '$schema': 'https://vega.github.io/schema/vega/v4.2.json',
    'width': 300,
    'height': 300,
    'data': [
        {'name': 'values',
         'values': [{'PC1': 0.0, 'PC2': 0.0}, {'PC1': 0.1, 'PC2': 0.2}, {'PC1': -0.5, 'PC2': -0.1}]},
    ],
    'scales': [
        {'name': 'xScale',
         'domain': {'data': 'values', 'field': 'PC1'},
         'range': 'width'},
        {'name': 'yScale',
         'domain': {'data': 'values', 'field': 'PC2'},
         'range': 'height'}],
    'axes': [
        {'scale': 'xScale', 'orient': 'bottom', 'title': 'PC1'},
        {'scale': 'yScale', 'orient': 'left', 'title': 'PC2'}],
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
                    "options": ["id1", "id2", "id3"],
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
                     'x': {'scale': 'xScale', 'field': 'PC1'},
                     'y': {'scale': 'yScale', 'field': 'PC2'}},
                 'update': {
                     'fill': [
                         {
                             'test': "datum.euclid_dist / datum.max_distance <= (1 - conservationLevel / 100)",
                             'value': '#3182bd'
                         },
                         {'value': 'black'}
                     ],
                     'opacity': [
                         {
                             'test': "datum[sequenceID] == null && hideMissingPositions",
                             'value': 0.0
                         },
                         {'value': 0.8}],
                     'tooltip': {
                         'signal': f"{{'title': 'position ' + datum['id'], "
                                   f"'PC1': datum['PC1'], "
                                   f"'PC2': datum['PC2']}}"}
                 }}}]}
