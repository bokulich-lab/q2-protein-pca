# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from skbio import OrdinationResults
from sklearn.decomposition import PCA


def _pca(ranks_df: pd.DataFrame,
         n_components: int = None) -> (OrdinationResults, OrdinationResults):
    # perform PCA
    pca_result = PCA(n_components=n_components)
    pca_result.fit(ranks_df)

    # transform ranks
    ranks_transformed = pd.DataFrame(pca_result.transform(ranks_df))
    ranks_transformed.index = ranks_df.index

    components_loadings = pd.DataFrame(-1 * pca_result.components_.T *
                                       np.sqrt(pca_result.explained_variance_))
    components_loadings.index = ranks_df.columns
    eigenvalues = pd.Series(pca_result.explained_variance_)

    ores_scores = OrdinationResults(
        short_method_name="PCA",
        long_method_name="Principal Components Analysis",
        eigvals=eigenvalues,
        samples=ranks_transformed,
        features=None,
        biplot_scores=None,
        proportion_explained=pd.Series(pca_result.explained_variance_ratio_))

    ores_loadings = OrdinationResults(
        short_method_name="PCA",
        long_method_name="Principal Components Analysis",
        eigvals=eigenvalues,
        samples=components_loadings,
        features=None,
        biplot_scores=None,
        proportion_explained=pd.Series(pca_result.explained_variance_ratio_))

    return ores_scores, ores_loadings


def pca(ranks: pd.DataFrame,
        n_components: int = None) -> (OrdinationResults, OrdinationResults):
    return _pca(ranks, n_components)
