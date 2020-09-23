from argparse import ArgumentParser

# Input arguments to the ProtAlign


def create_argument_parser():
    arg_parser = ArgumentParser(
        description='Fill out this description.')
    arg_parser.add_argument(
        '--verbosity',
        '-v',
        type=int,
        default=2,
        choices=range(1, 4),
        help='Verbosity level (1-3, default 2)')
    arg_parser.add_argument(
        '--main-location',
        type=str,
        help='Some main location.')
    arg_parser.add_argument(
        '--results-location',
        type=str,
        help='Results location.')
    arg_parser.add_argument(
        '--input-fasta-file',
        type=str,
        help='Name of the fasta file containing protein sequences.')
    arg_parser.add_argument(
        '--aln-slice-positions',
        type=str,
        help='Txt file with slicing positions for the alignment.')
    arg_parser.add_argument(
        '--clustalo-location',
        type=str,
        help='Location of the ClustalO executable.')
    arg_parser.add_argument(
        '--clustalo-overwrite',
        action="store_true",
        help='Should the files be overwritten.')
    arg_parser.add_argument(
        '--seq-start',
        type=int,
        default=0,
        help='Sequence start.')
    arg_parser.add_argument(
        '--seq-end',
        type=int,
        default=-1,
        help='Sequence end.')
    arg_parser.add_argument(
        '--expl-var-components',
        type=int,
        default=8,
        help='Number of PCA components to show.')
    arg_parser.add_argument(
        '--frac-pca-components',
        type=float,
        default=0.1,
        help='Fraction of PCA components to retain.')
    arg_parser.add_argument(
        '--components-to-project',
        type=int,
        default=2,
        help='Number of PCA components to project.')

    return arg_parser