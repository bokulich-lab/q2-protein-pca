import sys

from protein_pca.arguments import create_argument_parser
from protein_pca.alignment import ProtAlign
from protein_pca.helpers import set_logger
from protein_pca.pca import ProtPCA
from protein_pca.ranking import ProtRank


def run(args):
    logger = set_logger(args.verbosity, log_to_stdout=True)
    prot_align = ProtAlign(
        args.main_location,
        args.results_location,
        args.input_fasta_file,
        args.aln_slice_positions,
        args.clustalo_location,
        args.clustalo_overwrite,
        args.seq_start,
        args.seq_end,
    )
    prot_align.run()

    prot_rank = ProtRank(
        args.main_location, args.results_location, prot_align_obj=prot_align
    )
    prot_rank.run()

    prot_pca = ProtPCA(
        args.frac_pca_components,
        args.components_to_project,
        args.expl_var_components,
        prot_rank_obj=prot_rank,
    )
    prot_pca.run()


if __name__ == "__main__":
    arg_parser = create_argument_parser()
    args = arg_parser.parse_args(sys.argv[1:])
    run(args)
