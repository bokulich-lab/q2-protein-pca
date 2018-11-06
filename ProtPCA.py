import sys

from ProtPCA.arguments import create_argument_parser
from ProtPCA.alignment import ProtAlign
from ProtPCA.helpers import set_logger
from ProtPCA.pca import ProtPCA
from ProtPCA.ranking import ProtRank


def run(args):
    logger = set_logger(args.verbosity, log_to_stdout=True)
    prot_align = ProtAlign(args.config_file)
    prot_align.run()

    prot_rank = ProtRank(args.config_file, prot_align_obj=prot_align)
    prot_rank.run()

    prot_pca = ProtPCA(args.config_file, prot_rank_obj=prot_rank)
    prot_pca.run()

if __name__ == '__main__':
    arg_parser = create_argument_parser()
    args = arg_parser.parse_args(sys.argv[1:])
    run(args)