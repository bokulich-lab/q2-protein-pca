import logging
from datetime import datetime, timedelta

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from protein_pca import exceptions
from protein_pca.ranking import ProtRank

# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from pca_utils import eucl_dist, positions2color

logger = logging.getLogger(__name__)


def _adjust_ticks(axes):
    for ax in axes:
        for tick_x in ax.xaxis.get_major_ticks():
            tick_x.label.set_fontsize(14)
        for tick_xl in ax.xaxis.get_ticklabels():
            tick_xl.set_fontname('Arial')
        for tick_y in ax.yaxis.get_major_ticks():
            tick_y.label.set_fontsize(14)
        for tick_yl in ax.yaxis.get_ticklabels():
            tick_yl.set_fontname('Arial')


class ProtPCA:
    def __init__(self, frac_pca_components, components_to_project, expl_var_components, prot_rank_obj=None):
        self.timestamp_start = datetime.now()
        self.timestamp_end = None
        self.completion_time = None
        if prot_rank_obj is None:
            raise NotImplemented("You need to provide a valid ranked protein alignment.")
            # prot_rank_obj = ProtRank(config_file)
            # logger.warning('Ranked alignment object was not provided. Will attempt to read-in ranks from file.')
        else:
            self.main_loc = prot_rank_obj.main_loc
            self.all_seq_sliced = prot_rank_obj.all_seq_sliced
            self.summary_file = prot_rank_obj.summary_file
            self.clustal_out = prot_rank_obj.clustal_out
            self.clustal_out_fasta = prot_rank_obj.clustal_out_fasta
            self.log_file = prot_rank_obj.log_file
            self.alignment_length = prot_rank_obj.alignment_length
            self.alignment_df = prot_rank_obj.alignment_df
            self.alignment_df_copy = prot_rank_obj.alignment_df_copy
            self.occurrence_table = prot_rank_obj.occurrence_table
            self.ranking_table_out = prot_rank_obj.ranking_table_out
            if not prot_rank_obj.ranked:
                logger.error('Provided alignment is not ranked. Will try to read ranks from file.')
                self.copy_from_rank = False
                self.ranking_table = pd.read_csv(self.ranking_table_out, sep=',', header=0)
                if self.ranking_table.shape != [0,0]:
                    logger.info('Alignment ranks were successfully read from {}'.format(self.ranking_table_out))
                else:
                    logger.error('Rank import failed - empty file. Please check.')
                    raise exceptions.RankImportError('Rank import failed.')
            else:
                self.copy_from_rank = True
                self.ranking_table = prot_rank_obj.ranking_table
        self.components_to_keep = int(frac_pca_components*self.ranking_table.shape[1])
        self.components_to_project = components_to_project
        self.projections = dict()
        self.components = dict()
        self.pca = None
        self.fitted = False
        self.pca_transformed_data = None
        self.explained_variance = None
        self.explained_variance_ratio = None
        self.explained_variance_components = expl_var_components
        self.loadings = None

    def fit_pca(self):
        logger.info(f'PCA fitting... (restricted to {self.components_to_keep} components).')
        self.pca = PCA(n_components=self.components_to_keep)  # only keep the amount of principal components equal to 10% of all positions
        self.pca.fit(self.ranking_table)
        self.fitted = True
        self.pca_transformed_data = self.pca.transform(self.ranking_table)
        self.explained_variance = self.pca.explained_variance_
        self.explained_variance_ratio = self.pca.explained_variance_ratio_
        self.loadings = -1 * self.pca.components_.T * np.sqrt(self.explained_variance)
        self.timestamp_end = datetime.now()
        self.completion_time = (self.timestamp_end - self.timestamp_start) / timedelta(minutes=1)
        logger.info('Fitting completed in {:0.1f} minutes'.format(self.completion_time))

    def calculate_projections(self):
        logger.info(f'Calculating projections for {self.components_to_project} '
                    f'selected components (out of {self.components_to_keep} total).')
        if not self.fitted:
            logger.error('PCA was not fitted. Please check')
            raise exceptions.PCANotFitted()
        for i in range(self.components_to_project):
            self.projections[f'projection{i+1}'] = self.pca.components_[i]
            self.components[f'component{i+1}'] = self.pca_transformed_data[:, i]

    def plot_explained_variance(self, figsize=(8,5)):
        expl_variance_fig = plt.figure(figsize=figsize)
        variance_ax = expl_variance_fig.add_subplot(1, 1, 1)
        cumul_ax = expl_variance_fig.add_subplot(1, 1, 1)

        xvals = range(1, self.explained_variance_components + 1)
        yvals = self.explained_variance_ratio[:self.explained_variance_components]
        yvals_cum = np.cumsum(yvals)

        variance_ax.bar(xvals, yvals)
        cumul_ax.plot(xvals, yvals_cum, color='r', label='cumulative')

        cumul_ax.set_xlabel("Component", fontsize=16, fontname="Arial")
        cumul_ax.set_ylabel("Fraction explained variance", fontsize=16, fontname="Arial")
        cumul_ax.legend(loc="best", fontsize=14)
        cumul_ax.set_xticks(np.arange(1, self.explained_variance_components + 1))

        _adjust_ticks((cumul_ax))

        plt.tight_layout()
        return expl_variance_fig

    def plot_pca_scores(self, figsize=(8, 7)):
        scores_fig = plt.figure(figsize=figsize)

        scores_ax1 = scores_fig.add_subplot(2, 2, 1)  # PC1 vs PC2
        scores_ax2 = scores_fig.add_subplot(2, 2, 2)  # PC1 vs PC3
        scores_ax3 = scores_fig.add_subplot(2, 2, 3)  # PC2 vs PC3

        pc1_data = self.pca_transformed_data[:,0]
        pc2_data = self.pca_transformed_data[:,1]
        pc3_data = self.pca_transformed_data[:,2]

        scores_ax1.scatter(pc1_data, pc2_data, c='r', s=30)
        scores_ax2.scatter(pc1_data, pc3_data, c='b', s=30)
        scores_ax3.scatter(pc2_data, pc3_data, c='g', s=30)

        scores_ax1.set_xlabel("PC1", fontsize=16, fontname="Arial")
        scores_ax1.set_ylabel("PC2", fontsize=16, fontname="Arial")
        scores_ax2.set_xlabel("PC1", fontsize=16, fontname="Arial")
        scores_ax2.set_ylabel("PC3", fontsize=16, fontname="Arial")
        scores_ax3.set_xlabel("PC2", fontsize=16, fontname="Arial")
        scores_ax3.set_ylabel("PC3", fontsize=16, fontname="Arial")

        _adjust_ticks((scores_ax1, scores_ax2, scores_ax3))

        #plt.suptitle('PCA scores', fontsize=20)
        plt.tight_layout(pad=1, h_pad=1.5, w_pad=1.5)

        return scores_fig

    def run(self):
        self.fit_pca()
        self.calculate_projections()
        some_fig = self.plot_pca_scores()
        some_fig.savefig('/Users/misialq/Desktop/testfigure.pdf', dpi=300)
