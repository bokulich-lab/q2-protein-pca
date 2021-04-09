# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import skbio
import skbio.io
from q2_types.feature_data._transformer import AlignedProteinIterator

from q2_types.feature_data import ProteinFASTAFormat, AlignedProteinFASTAFormat
from q2_alignment._mafft import run_command


def _mafft(sequences_fp, alignment_fp, n_threads, parttree):
    # Save original sequence IDs since long ids (~250 chars) can be truncated
    # by mafft. We'll replace the IDs in the aligned sequences file output by
    # mafft with the originals.
    #
    # https://github.com/qiime2/q2-alignment/issues/37
    aligned_seq_ids = {}
    unaligned_seq_ids = {}

    # if alignment_fp is not None:
    #     for seq in skbio.io.read(alignment_fp, format='fasta',
    #                              constructor=skbio.Protein):
    #         id_ = seq.metadata['id']
    #         if id_ in aligned_seq_ids:
    #             raise ValueError(
    #                 "A sequence ID is duplicated in the aligned sequences: "
    #                 "%r" % id_)
    #         else:
    #             aligned_seq_ids[id_] = True

    for seq in skbio.io.read(sequences_fp, format='fasta',
                             constructor=skbio.Protein):
        id_ = seq.metadata['id']
        if id_ in unaligned_seq_ids:
            raise ValueError(
                "A sequence ID is duplicated in the unaligned sequences: "
                "%r" % id_)
        elif id_ in aligned_seq_ids:
            raise ValueError(
                "A sequence ID is present in both the aligned and unaligned "
                "sequences: %r" % id_)
        else:
            unaligned_seq_ids[id_] = True

    result = AlignedProteinFASTAFormat()
    result_fp = str(result)
    ids = {**aligned_seq_ids, **unaligned_seq_ids}

    # mafft will fail if the number of sequences is larger than 1 million.
    # mafft requires using parttree which is an algorithm to build an
    # approximate tree from a large number of unaligned sequences.
    # By catching the error below if a user has not used parttree flag, we are
    # eliminating the need for the mafft error to be shown to the user which
    # can be confusing and intimidating.

    if not parttree and len(ids) > 1000000:
        raise ValueError(
            "The number of sequences in your feature table is larger than "
            "1 million, please use the parttree parameter")

    # mafft's signal for utilizing all cores is -1. We want to our users
    # to enter auto for using all cores. This is to prevent any confusion and
    # to keep the UX consisent.
    if n_threads == 'auto':
        n_threads = -1

    # `--inputorder` must be turned on because we need the input and output in
    # the same sequence order to replace the IDs below. This is mafft's default
    # behavior but we pass the flag in case that changes in the future.
    cmd = ["mafft", "--preservecase", "--inputorder",
           "--thread", str(n_threads)]

    if parttree:
        cmd += ['--parttree']

    if alignment_fp is not None:
        cmd += ['--add', sequences_fp, alignment_fp]
    else:
        cmd += [sequences_fp]

    run_command(cmd, result_fp)

    # Read output alignment into memory, reassign original sequence IDs, and
    # write alignment back to disk.
    msa = skbio.TabularMSA.read(result_fp, format='fasta',
                                constructor=skbio.Protein)
    # Using `assert` because mafft would have had to add or drop sequences
    # while aligning, which would be a bug on mafft's end. This is just a
    # sanity check and is not expected to trigger in practice.
    assert len(ids) == len(msa)
    for id, seq in zip(ids, msa):
        seq.metadata['id'] = id

    # Turning off roundtripping options to speed up writing. We can safely turn
    # these options off because we know the sequence IDs are rountrip-safe
    # since we read them from a FASTA file above.
    #
    # http://scikit-bio.org/docs/latest/generated/
    #     skbio.io.format.fasta.html#writer-specific-parameters
    msa.write(result_fp, id_whitespace_replacement=None,
              description_newline_replacement=None)
    return result


def _map_positions(aligned_sequences: AlignedProteinIterator) -> pd.DataFrame:
    mapping_df = pd.DataFrame()
    for aln_seq_record in aligned_sequences:
        id_aln = aln_seq_record.metadata['id']
        seq_aln = pd.Series(aln_seq_record.values.astype('str'), name=id_aln)

        seq_aln_degapped = seq_aln[seq_aln != "-"]
        original_positions = pd.Series([int(x) for x in range(
            len(seq_aln_degapped.index))], index=seq_aln_degapped.index)
        seq_aln[seq_aln != "-"] = original_positions
        seq_aln.replace(to_replace="-", value=np.nan, inplace=True)
        seq_aln = seq_aln.astype('Int64')

        mapping_df = pd.concat([mapping_df, seq_aln], axis=1)

    mapping_df.index.name = 'Alignment position'
    mapping_df.index = mapping_df.index.astype('int64')

    return mapping_df


def mafft(sequences: ProteinFASTAFormat,
          n_threads: int = 1,
          parttree: bool = False) -> AlignedProteinFASTAFormat:
    sequences_fp = str(sequences)
    return _mafft(sequences_fp, None, n_threads, parttree)


def map_positions(
        aligned_sequences: AlignedProteinIterator) -> pd.DataFrame:
    return _map_positions(aligned_sequences)
