#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_protein_pca._format import (
    PositionMappingFormat, PositionMappingDirectoryFormat,
    RankedProteinAlignmentDirectoryFormat, RankedProteinAlignmentFormat)
from q2_protein_pca._type import PositionMapping, RankedProteinAlignment
from q2_types.feature_data._type import (
    ProteinSequence, AlignedProteinSequence, FeatureData)
from q2_types.ordination import PCoAResults
from qiime2.plugin import (Str, Plugin, Choices, Bool, Citations, Int, Range)

import q2_protein_pca

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
    input_descriptions={'sequences': 'Protein sequences to be aligned.'},
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use `auto` to automatically use '
                     'all available cores)',
        'parttree': 'This flag is required if the number of sequences being '
                    'aligned are larger than 1000000. Disabled by default'},
    output_descriptions={'alignment': 'Aligned protein sequences.'},
    name='De novo multiple protein sequence alignment with MAFFT',
    description=(
        "Perform de novo multiple protein sequence alignment using MAFFT."),
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
    description=(
        "Perform protein alignment ranking based on amino acid "
        "occurrence frequency."),
    citations=[citations['Wang2014']]
)

plugin.methods.register_function(
    function=q2_protein_pca.pca,
    inputs={'ranks': FeatureData[RankedProteinAlignment]},
    parameters={'n_components': Int % Range(1, None)},
    outputs=[('pca_scores', PCoAResults), ('pca_loadings', PCoAResults)],
    input_descriptions={'ranks': 'Ranked protein alignment.'},
    parameter_descriptions={
        'n_components': 'The number of principal components to retain.'},
    output_descriptions={
        'pca_scores': 'PCA scores.',
        'pca_loadings': 'PCA loadings.'},
    name='Principal Component Analysis of ranked protein alignment',
    description=(
        "Perform PCA on protein alignment ranked according to amino acid "
        "occurence frequencies."),
    citations=[citations['Wang2014']]
)

plugin.methods.register_function(
    function=q2_protein_pca.map_positions,
    inputs={'aligned_sequences': FeatureData[AlignedProteinSequence]},
    parameters={},
    outputs=[('mapped_positions', FeatureData[PositionMapping])],
    input_descriptions={'aligned_sequences': 'Aligned protein sequences.'},
    parameter_descriptions={},
    output_descriptions={
        'mapped_positions': 'Amino acid positions mapping between raw '
                            '(unaligned) and aligned sequences.'},
    name='Amino acid position mapping',
    description=(
        "Find mapping of amino acid positions between an aligned protein "
        "sequence and its unaligned counterpart."),
)

plugin.visualizers.register_function(
    function=q2_protein_pca.plot_loadings,
    inputs={'pca_loadings': PCoAResults,
            'positions_mapping': FeatureData[PositionMapping]},
    parameters={'pdb_id': Str, 'nterm_offset': Int % Range(0, None)},
    input_descriptions={'pca_loadings': 'PCA loadings.',
                        'positions_mapping': 'Amino acid positions mapping.'},
    parameter_descriptions={'pdb_id': 'PDB ID of the protein structure to '
                                      'display conserved positions on.',
                            'nterm_offset': 'Number of the amino acids that'
                                            'are missing from the N-terminus'
                                            'of the PDB structure. Defaults'
                                            'to 1.'},
    name='PCA loadings plot',
    description=(
        'Visualise principal component loadings to find which positions '
        'within the alignment/protein sequence contribute most/least to'
        'the observed variance. Ultimately, find which positions are most/'
        'least conserved within a protein sequence.')
)

# Registrations
plugin.register_formats(PositionMappingFormat, PositionMappingDirectoryFormat)
plugin.register_formats(
    RankedProteinAlignmentFormat,
    RankedProteinAlignmentDirectoryFormat)

plugin.register_semantic_types(PositionMapping)
plugin.register_semantic_types(RankedProteinAlignment)

plugin.register_semantic_type_to_format(
    FeatureData[PositionMapping],
    artifact_format=PositionMappingDirectoryFormat)
plugin.register_semantic_type_to_format(
    FeatureData[RankedProteinAlignment],
    artifact_format=RankedProteinAlignmentDirectoryFormat)

importlib.import_module('q2_protein_pca._transformer')
