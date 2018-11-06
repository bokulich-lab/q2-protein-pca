import os
import pprint

import pandas as pd
import numpy as np
import logging
from datetime import datetime, timedelta

from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
#from pca_utils import eucl_dist, positions2color

from ProtPCA import exceptions
from ProtPCA.ranking import ProtRank

logger = logging.getLogger()


class ProtPCA:
    def __init__(self, config_file, prot_rank_obj=None):
        self.timestamp_start = datetime.now()
        self.timestamp_end = None
        self.completion_time = None
        if prot_rank_obj is None:
            prot_rank_obj = ProtRank(config_file)
            logger.warning('Ranked alignment object was not provided. Will attempt to read-in ranks from file.')
        else:
            self.config = prot_rank_obj.config
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
        self.components_to_keep = int(self.config['analysis']['frac_pca_components']*self.ranking_table.shape[1])
        self.components_to_project = self.config['analysis']['components_to_project']
        self.projections = dict()
        self.components = dict()
        self.pca = None
        self.fitted = False
        self.pca_transformed_data = None

    def fit_pca(self):
        logger.info(f'PCA fitting... (restricted to {self.components_to_keep} components).')
        self.pca = PCA(n_components=self.components_to_keep)  # only keep the amount of principal components equal to 10% of all positions
        self.pca.fit(self.ranking_table)
        self.fitted = True
        self.pca_transformed_data = self.pca.transform(self.ranking_table)
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

    def run(self):
        self.fit_pca()
        self.calculate_projections()
