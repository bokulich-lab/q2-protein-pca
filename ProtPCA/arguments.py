from argparse import ArgumentParser

# Input arguments to the ProtAlign


def create_argument_parser():
    arg_parser = ArgumentParser(
        description='Utility that listens to a Redis channel for audio, '
                    'and dumps the data as separate audio files to a '
                    'specified location.')
    arg_parser.add_argument(
        '-cf',
        '--config-file',
        type=str,
        help='Path to the config file')
    arg_parser.add_argument(
        '--verbosity',
        '-v',
        type=int,
        default=2,
        choices=range(1, 4),
        help='Verbosity level (1-3, default 2)')

    return arg_parser