import os
import pprint
import logging
from tqdm import tqdm
from datetime import datetime, timedelta

import configobj
import pandas as pd
import pkg_resources
import validate
from Bio import SeqIO

from ProtPCA import exceptions
from ProtPCA.helpers import split_string

logger = logging.getLogger()


class ProtRank:
    def __init__(self, config_file, prot_align_obj=None):
        self.timestamp_start = datetime.now()
        self.timestamp_end = None
        self.completion_time = None
        self.ranked = False
        if prot_align_obj is None:
            # load and validate the config
            self.cwd = os.path.dirname(os.path.realpath(__file__))
            validator = validate.Validator()
            if config_file is None:
                self.config = {}
            self.config = configobj.ConfigObj(
                config_file,
                configspec=pkg_resources.resource_filename(
                    __name__,
                    '../config/ProtPCA.spec'))
            self.validated = self.config.validate(validator)
            if isinstance(self.validated, dict) or not self.validated:
                logger.error(pprint.pformat(self.validated))
                raise exceptions.ConfigNotValid('Configuration file contains error')

            self.main_loc = self.config['data']['inputs']['main_loc']
            if not os.path.isdir(os.path.join(self.main_loc, self.config['data']['outputs']['results_folder'])):
                logger.warning('Path {} does not exist. Creating...'.format(
                    os.path.join(self.main_loc, self.config['data']['outputs']['results_folder'])
                ))
                os.mkdir(os.path.join(self.main_loc, self.config['data']['outputs']['results_folder']))
            self.all_seq_sliced = os.path.join(
                self.main_loc,
                self.config['data']['outputs']['results_folder'],
                '{}_all_sequences_sliced.fasta'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
            )
            self.summary_file = os.path.join(
                self.main_loc,
                self.config['data']['outputs']['results_folder'],
                '{}_final_summary.txt'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
            )
            self.clustal_out = os.path.join(
                self.main_loc,
                self.config['data']['outputs']['results_folder'],
                '{}_sliced_alignment.aln'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
            )
            self.clustal_out_fasta = os.path.join(
                self.main_loc,
                self.config['data']['outputs']['results_folder'],
                '{}_sliced_alignment.fasta'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
            )
            self.log_file = os.path.join(
                self.main_loc,
                self.config['data']['outputs']['results_folder'],
                '{}_log.txt'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
            )
            self.copy_from_aln = False
        else:
            self.config = prot_align_obj.config
            self.main_loc = prot_align_obj.main_loc
            self.all_seq_sliced = prot_align_obj.all_seq_sliced
            self.summary_file = prot_align_obj.summary_file
            self.clustal_out = prot_align_obj.clustal_out
            self.clustal_out_fasta = prot_align_obj.clustal_out_fasta
            self.log_file = prot_align_obj.log_file
            self.copy_from_aln = True
        self.alignment_length = 0
        self.alignment_df = None
        self.alignment_df_copy = None
        self.occurrence_table = None
        self.ranking_table = None
        self.ranking_table_out = os.path.join(
                self.main_loc,
                self.config['data']['outputs']['results_folder'],
                '{}_ranking_table.csv'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
            )

    def read_alignment(self, aln_file, aln_type='fasta'):
        logger.info('Reading in the alignment file...')
        master_list = [split_string(str(seq_record.seq)) for seq_record in SeqIO.parse(aln_file, aln_type)]
        self.alignment_length = len(master_list[0])
        col_names = [f'pos{x}' for x in range(self.alignment_length)]
        self.alignment_df = pd.DataFrame(master_list, columns=col_names)
        self.alignment_df_copy = self.alignment_df.copy(deep=True)
        logger.info(f'Alignment length: {self.alignment_length}, total sequence count: {self.alignment_df.shape[0]}')

    @staticmethod
    def get_occurrences(df):
        return df.apply(pd.value_counts).fillna(0).astype('int64')

    @staticmethod
    def rank_alphabet(pd_series):
        # remove gaps beofore ranking
        pd_series_1 = pd_series.loc[[x for x in pd_series.index if x is not '-']]
        # sort AAs according to occurrence (descending)
        temp_df = pd.DataFrame({'idx_col': pd_series_1.index, pd_series_1.name: pd_series_1})
        test_rank_x = temp_df.sort_values(by=[pd_series_1.name, 'idx_col'], ascending=[False, False])
        test_rank = pd.Series(data=test_rank_x[pd_series_1.name].values, index=test_rank_x['idx_col'].values)

        # group by count of occurrences
        all_occurrences = test_rank.groupby(test_rank).agg('count')
        # list how many AAs with given count (no of occur, no of aas) in a count descending order
        list_of_rep_aas = sorted(list(all_occurrences.iteritems()), key=lambda x: x[0],
                                 reverse=True)
        # generates a list of tuples containing the amino acid and its occurrence frequency, e.g.:
        # [('A', 10), ('C', 10), ('S', 0), ('R', 0), ('B', 0), ('D', 0) ... ]
        all_items = list(
            test_rank.iteritems())
        all_items_ranked = []
        all_ranks_by_sum = 0
        for i, e in enumerate(list_of_rep_aas):
            all_ranks_by_sum += e[1] if e[0] != 0 else False

        cnt1 = 0
        tot_aas = all_ranks_by_sum
        for x, rep_aa in enumerate(list(list_of_rep_aas)):
            cnt2 = 0
            rep_aa_count, rep_aa_aanum = rep_aa
            if rep_aa_count != 0.0:
                list_of_ranks = []
                for i in range(rep_aa_aanum):
                    rank = tot_aas - rep_aa_aanum + cnt2 + 1
                    list_of_ranks.append(rank)
                    cnt2 += 1
                tot_aas = tot_aas - rep_aa_aanum

                for j, j_rank in enumerate(list_of_ranks):
                    aas_with_rank = (all_items[cnt1][0], j_rank)
                    all_items_ranked.append(aas_with_rank)
                    cnt1 += 1
        ranks_dict = dict(all_items_ranked)

        # replace gap (-) counts with 0
        ranks_dict.update({'-': 0})

        return ranks_dict

    def rank_all_columns(self, aln_df):
        logger.info('Ranking alignment positions...')
        # calculate AA occurrences
        occur_freq = self.get_occurrences(aln_df)
        self.occurrence_table = occur_freq
        # calculate AA ranks
        tqdm.pandas()
        aa_ranks = occur_freq.progress_apply(self.rank_alphabet, axis=0)
        # replace original AAs with their ranks
        aln_df_ranked = aln_df.replace(aa_ranks)
        self.ranked = True
        logger.info('Ranking completed.')
        return aln_df_ranked

    def save_ranks_to_file(self):
        if self.ranking_table is not None:
            self.ranking_table.to_csv(self.ranking_table_out, sep=',', index=False, header=True)

    def run(self):
        logger.info('Initiating ranking...')
        self.read_alignment(self.clustal_out_fasta, 'fasta')
        self.ranking_table = self.rank_all_columns(self.alignment_df)
        self.save_ranks_to_file()
        self.timestamp_end = datetime.now()
        self.completion_time = (self.timestamp_end - self.timestamp_start) / timedelta(minutes=1)
        logger.info('Alignment ranking completed in {:0.1f} minutes'.format(self.completion_time))


