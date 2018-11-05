import os
import pprint
import logging
from datetime import datetime, timedelta

import configobj
import validate
import pandas as pd
import pkg_resources
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

from . import exceptions

logger = logging.getLogger()


class ProtAlign:
    def __init__(self, config_file):
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

        self.timestamp_start = datetime.now()
        self.timestamp_end = None
        self.completion_time = None
        self.main_loc = self.config['data']['inputs']['main_loc']
        self.all_seq_fasta = os.path.join(self.main_loc, self.config['data']['inputs']['input_fasta'])
        self.positions_file = os.path.join(self.main_loc, self.config['data']['inputs']['aln_slice_positions'])
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
        self.clustalo_exe = r"{}".format(self.config['data']['inputs']['clustalo_location'])
        self.log_file = os.path.join(
            self.main_loc,
            self.config['data']['outputs']['results_folder'],
            '{}_log.txt'.format(self.timestamp_start.strftime("%Y%m%d-%H%M%S"))
        )
        self.ClustalForceOverwrite = self.config['data']['inputs']['clustalo_overwrite']

        self.record_iterator = SeqIO.parse(self.all_seq_fasta, "fasta")
        self.ids, self.seqs = list(), list()
        for record in self.record_iterator:
            self.ids.append(record.id)
            self.seqs.append(record.seq)
        self.seq_lens = [len(seq) for seq in self.seqs]
        self.sequences_df = pd.DataFrame({'ID': self.ids, 'Sequence': self.seqs, 'Length': self.seq_lens})
        self.max_seq = max(self.seq_lens)
        self.min_seq = min(self.seq_lens)
        logger.info("The longest sequence is {0} AA long. The shortest sequence is {1} AA long."
                    .format(self.max_seq, self.min_seq))

        self.sequence_start = self.config['data']['inputs']['sequence_start']
        self.sequence_end = self.max_seq if self.config['data']['inputs']['sequence_end'] == -1 else \
            self.config['data']['inputs']['sequence_end']

        self.positions = None

    def read_aln_positions(self):
        """Reads positions file and prepares a list of slices for every sequence"""
        try:
            with open(self.positions_file, 'r') as f:
                lines = f.readlines()
                self.positions = [tuple(x.strip().split(',')) for x in lines]
            if not self.positions:
                logger.warning("No positions could be found. Sequences will be sliced to {}-{}.".format(
                    self.sequence_start, self.sequence_end
                ))
                self.positions = [(self.sequence_start, self.sequence_end) for id in self.ids]
        except IOError:
            logger.warning("'Positions' file ({}) not found. Sequences will be sliced to {}-{}.".format(
                self.positions_file, self.sequence_start, self.sequence_end
            ))
            self.positions = [(self.sequence_start, self.sequence_end) for id in self.ids]
        finally:
            logger.info("List of positions has {0} entries.".format(len(self.positions)))
            if len(self.ids) == len(self.seqs) == len(self.positions):
                logger.info("Found and successfully processed {0} sequences.".format(len(self.ids)))
                self.sequences_df['Slice_start'] = [x[0] for x in self.positions]
                self.sequences_df['Slice_end'] = [x[1] for x in self.positions]
                self.sequences_df['Slice_length'] = self.sequences_df['Slice_end'] - self.sequences_df['Slice_start']
            else:
                logger.error("Mismatch found. {0} IDs {1} sequences and {2} positions were found."
                             .format(len(self.ids), len(self.seqs), len(self.positions)))
                raise exceptions.SequenceMismatch('Mismatch between number of sequences, ids and positions.')

    @staticmethod
    def slice_sequence(sequence, start, end):
        try:
            pos0, pos1 = int(start.strip()) - 1, int(end.strip())
        except AttributeError:
            pos0, pos1 = int(start) - 1, int(end)
        return sequence[pos0:pos1] if pos1 <= len(sequence) else sequence

    def generate_sequence_slices(self):
        self.sequences_df['Sliced_sequence'] = self.sequences_df.apply(
            lambda row: self.slice_sequence(row['Sequence'], row['Slice_start'], row['Slice_end']), axis=1)

    def write_slices_to_fasta(self):
        """Writes sliced sequences to a new multi-line FASTA file"""
        with open(self.all_seq_sliced, "w") as f:
            for i, row in self.sequences_df.iterrows():
                f.write(f'>{row["ID"]}\n{row["Sliced_sequence"]}\n\n')

    def align_sequences(self):
        logger.info('Alignment is starting...')
        clustalo_cline = ClustalOmegaCommandline(self.clustalo_exe,
                                                 infile=self.all_seq_sliced,
                                                 outfile=self.clustal_out,
                                                 force=self.ClustalForceOverwrite,
                                                 verbose=True,
                                                 auto=True)
        logger.info(clustalo_cline)

        if os.path.isfile(self.clustalo_exe):
            stdout, stderr = clustalo_cline()
        else:
            logger.error('ClustalO executable could not be found.')
            raise exceptions.MissingClustalOExecutable('ClustalO executable could not be found.')
        logger.info("ClustalO alignment successfully completed.")

    @staticmethod
    def convert_aln_to_fasta(aln_file, fasta_output):
        """convert ALN alignment file into FASTA format"""
        with open(aln_file, 'r') as aln_f:
            aln_fasta_file = open(fasta_output, 'w')
            for line in aln_f.readlines():
                aln_fasta_file.write(line)
            aln_fasta_file.close()
        logger.info(f'File {aln_file} successfully converted to {fasta_output}.')

    def run(self):
        logger.info('Initiating alignments...')
        self.read_aln_positions()
        self.generate_sequence_slices()
        self.write_slices_to_fasta()
        self.align_sequences()
        self.convert_aln_to_fasta(self.clustal_out, self.clustal_out_fasta)
        self.timestamp_end = datetime.now()
        self.completion_time = (self.timestamp_end - self.timestamp_start) / timedelta(minutes=1)
        logger.info('Sequence alignment completed in {:0.1f} minutes'.format(self.completion_time))
