[data]
    [[inputs]]
        main_loc = string()
        input_fasta = string()
        aln_slice_positions = string(default='positions.txt')
        clustalo_location = string(default='~/clustalo-1.2.4')
        clustalo_overwrite = boolean(default=True)
        sequence_start = integer(min=0, default=0)
        sequence_end = integer(default=-1)
    [[outputs]]
        results_folder = string(default='analysis')