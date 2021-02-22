#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

setup(
    name='q2-protein-pca',
    version='2020.08',
    license='BSD-3-Clause',
    packages=find_packages(),
    install_requires=['biopython>=1.78', 'scikit-learn', 'pandas',
                      'numpy', 'matplotlib'],
    author="Michal Ziemski",
    author_email="ziemski.michal@gmail.com",
    description=("Dimensionality reduction of protein sequences using "
                 "PCA on ranked sequence alignments."),
    url="https://github.com/misialq/pca-protein-analysis",
    entry_points={
        'qiime2.plugins':
        ['q2-protein-pca=q2_protein_pca.plugin_setup:plugin']
    },
    package_data={
        'q2_protein_pca': [
            'data/index.html',
            'citations.bib',
            'assets/loadings/index.html',
            'assets/loadings/vega/css/*',
            'assets/loadings/vega/js/*',
            'assets/loadings/vega/licenses/*',
            'assets/loadings/pdb-litemol/*',
            'assets/loadings/pdb-litemol/license/*',
            'assets/loadings/assets/fonts/*'
        ],
        'q2_protein_pca.tests': ['data/*'],
    },
    zip_safe=False,
)
