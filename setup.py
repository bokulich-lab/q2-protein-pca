# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name='q2-protein-pca',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Michal Ziemski",
    author_email="ziemski.michal@gmail.com",
    description=("Dimensionality reduction of protein sequences using "
                 "PCA on ranked sequence alignments."),
    url="https://github.com/bokulich-lab/q2-protein-pca",
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
