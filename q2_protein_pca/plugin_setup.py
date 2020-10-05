#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.feature_data import FeatureData
from q2_types.feature_data._type import ProteinSequence, AlignedProteinSequence, RankedProteinAlignment
from q2_types.ordination import PCoAResults
from qiime2.plugin import (Str, Plugin, Metadata, Choices, Bool, Citations,
                           Int, MetadataColumn, Numeric, Range)

import q2_protein_pca
import importlib
from q2_types.sample_data import SampleData
from q2_types.distance_matrix import DistanceMatrix
from q2_types.tree import Phylogeny, Rooted

citations = Citations.load('citations.bib', package='q2_protein_pca')

plugin = Plugin(
    name='protein-pca',
    version=q2_protein_pca.__version__,
    website="https://github.com/misialq/pca-protein-analysis",
    package='q2_protein_pca',
    description=(
        'This QIIME 2 plugin supports methods for PCA analysis of ranked'
        ' protein alignments.'),
    short_description=(
        'Plugin for PCA analysis of protein sequences.'),
)

plugin.methods.register_function(
    function=q2_protein_pca.mafft,
    inputs={'sequences': FeatureData[ProteinSequence]},
    parameters={'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
                'parttree': Bool},
    outputs=[('alignment', FeatureData[AlignedProteinSequence])],
    input_descriptions={'sequences': 'The sequences to be aligned.'},
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use `auto` to automatically use '
                     'all available cores)',
        'parttree': 'This flag is required if the number of sequences being '
                    'aligned are larger than 1000000. Disabled by default'},
    output_descriptions={'alignment': 'The aligned protein sequences.'},
    name='De novo multiple protein sequence alignment with MAFFT',
    description=("Perform de novo multiple protein sequence alignment using MAFFT."),
    citations=[citations['katoh2013mafft']]
)

plugin.methods.register_function(
    function=q2_protein_pca.rank_alignment,
    inputs={'sequences': FeatureData[AlignedProteinSequence]},
    parameters={},
    outputs=[('ranked_alignment', FeatureData[RankedProteinAlignment])],
    input_descriptions={'sequences': 'Aligned protein sequences.'},
    parameter_descriptions={},
    output_descriptions={'ranked_alignment': 'Ranked protein alignment.'},
    name='Protein alignment ranking',
    description=("Perform protein alignment ranking based on amino acid "
                 "occurence frequency."),
    citations=[citations['Wang2014']]
)

plugin.methods.register_function(
    function=q2_protein_pca.pca,
    inputs={'ranks': FeatureData[RankedProteinAlignment]},
    parameters={'n_components': Int % Range(1, None)},
    outputs=[('pca_scores', PCoAResults), ('pca_loadings', PCoAResults)],
    input_descriptions={'ranks': 'Ranked protein alignment.'},
    parameter_descriptions={'n_components': 'The number of principal '
                                            'components to retain.'},
    output_descriptions={'pca_scores': 'PCA scores.', 'pca_loadings': 'PCA loadings.'},
    name='Principal Compnent Analysis of ranked protein alignment',
    description=("Perform PCA on protein alignment ranked according to "
                 "amino acid occurence frequencies."),
    citations=[citations['Wang2014']]
)

plugin.visualizers.register_function(
    function=q2_protein_pca.plot_loadings,
    inputs={'pca_loadings': PCoAResults},
    parameters={},
    input_descriptions={'pca_loadings': 'PCA loadings.'},
    parameter_descriptions={},
    name='Test visualizer',
    description=('Test visualizer description.')
)

# base_parameters = {
#     'metadata': Metadata,
#     'latitude': Str,
#     'longitude': Str,
#     'missing_data': Str
# }
#
# base_parameter_descriptions = {
#     'metadata': 'The sample metadata containing coordinate data.',
#     'latitude': 'Metadata column containing latitude in decimal degrees.',
#     'longitude': 'Metadata column containing longitude in decimal degrees.',
#     'missing_data': 'If "error" (default), will raise an error if any '
#                     'metadata columns are missing data. Set to "ignore" to '
#                     'silently drop rows (samples) that are missing data.'
# }


# plugin.visualizers.register_function(
#     function=draw_map,
#     inputs={},
#     parameters={**base_parameters,
#                 'column': Str,
#                 'color_palette': Str % Choices([
#                     'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired',
#                     'Accent', 'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c',
#                     'viridis', 'plasma', 'inferno', 'magma', 'terrain',
#                     'rainbow']),
#                 'discrete': Bool,
#                 'image': Str % Choices(
#                     ['StamenTerrain', 'OSM', 'GoogleTiles']),
#                 },
#     input_descriptions={},
#     parameter_descriptions={
#         **base_parameter_descriptions,
#         'column': ('Metadata column to use for coloring sample points. If '
#                    'none is supplied, will use alpha_diversity artifact for '
#                    'coloring.'),
#         'color_palette': (
#             'Color palette to use for coloring sample points on map.'),
#         'discrete': 'Plot continuous column data as discrete values.',
#         'image': 'Base map image to use for coordinate projection.'},
#     name='Plot sampling site geocoordinates on a map.',
#     description=('Plots sample data onto a map using sample geocoordinates. '
#                  'Sample points are colored by the column name "column", '
#                  'which may be categorical or numeric. Note that samples with '
#                  'missing values are silently dropped.'),
#     citations=[citations['Cartopy']]
# )
#
# plugin.methods.register_function(
#     function=geodesic_distance,
#     inputs={},
#     parameters=base_parameters,
#     outputs=[('distance_matrix', DistanceMatrix)],
#     input_descriptions={},
#     parameter_descriptions=base_parameter_descriptions,
#     name='Create a distance matrix from sample geocoordinates.',
#     description='Measure pairwise geodesic distances between coordinates. '
#                 'Output distances are reported in meters. '
#                 'Note that samples with missing values are silently dropped.',
#     citations=[citations['Karney2013']]
# )

# Registrations
# plugin.register_formats(CoordinatesFormat, CoordinatesDirectoryFormat)
#
# plugin.register_semantic_types(Coordinates)
#
# plugin.register_semantic_type_to_format(
#     SampleData[Coordinates],
#     artifact_format=CoordinatesDirectoryFormat)
#
# plugin.register_formats(QuadTreeFormat, QuadTreeDirectoryFormat)
#
# plugin.register_semantic_types(QuadTree)
#
# plugin.register_semantic_type_to_format(
#     SampleData[QuadTree],
#     artifact_format=QuadTreeDirectoryFormat)
# importlib.import_module('q2_coordinates._transformer')
