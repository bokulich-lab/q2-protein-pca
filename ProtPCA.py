import sys

from ProtPCA.arguments import create_argument_parser
from ProtPCA.alignment import ProtAlign
from ProtPCA.helpers import set_logger


def run(args):
    logger = set_logger(args.verbosity, log_to_stdout=True)
    prot_align = ProtAlign(args.config_file)
    prot_align.run()


if __name__ == '__main__':
    arg_parser = create_argument_parser()
    args = arg_parser.parse_args(sys.argv[1:])
    run(args)