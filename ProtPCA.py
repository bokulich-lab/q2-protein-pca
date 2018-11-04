from ProtPCA import alignment
import sys

from ProtPCA.arguments import create_argument_parser
from ProtPCA.alignment import ProtAlign


def run(args):
    prot_align = ProtAlign(args.config_file)
    prot_align.run()



if __name__ == '__main__':
    arg_parser = create_argument_parser()
    args = arg_parser.parse_args(argv[1:])
    run(args)