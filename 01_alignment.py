from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline


from datetime import datetime
from time import time
import os

### align many sequences using a list of specific positions
time0 = time()

#####################################-----------------VARIABLES----------------##################################
# indicate main subfolders where all the data will be saved:
subfolder = "trx_sequences/analysis" # -> location of the FASTA file containing all the sequences to be aligned
subfolder2 = "trx_sequences" # -> should contain all the analysis results at the end

# (in) FASTA file with all protein sequences
all_seq_fasta = subfolder2 + "/trx_sequences_all.fasta"

# (in) location of slicing positions (in case you do not want to align full-length sequences)
positions_file = subfolder2 + "/positions_txt.txt"

# (out) FASTA file for saving sequences after slicing
all_seq_sliced = subfolder + "/all_sequences_sliced_" + datetime.now().strftime(
    "%Y%m%d-%H%M%S") + ".fasta"

# (out) locations for alignment files (from clustal)
clustal_out = subfolder + "/sliced_alignment_" + datetime.now().strftime(
    "%Y%m%d-%H%M%S") + ".aln"
clustal_out_fasta = subfolder + "/sliced_alignment_" + datetime.now().strftime(
    "%Y%m%d-%H%M%S") + ".fasta"

# (out) location for a text file containing summary of all sequences
summary_file = subfolder + "/final_summary_" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".txt"

# clustal exe location
clustalo_exe = r"/clustalo-1.2.4"

# (out) location of a log file
log_file = subfolder + "/log_" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".txt"

# force file overwrite by clustal
forceO = True

###############################################################################################################

log = open(log_file, "w")
log.write("Log file generated on: {0}\n\n".format(str(datetime.now())))

### generate an iterator for all sequences from a fasta file
record_iterator = SeqIO.parse(all_seq_fasta, "fasta")

### loop through all the entries in the fasta file and extract ids & sequences
ids = []
seqs = []
for record in record_iterator:
    ids.append(record.id)
    seqs.append(record.seq)

### find the longest sequence
seq_lens = [len(seq) for seq in seqs]
max_seq = max(seq_lens)
min_seq = min(seq_lens)
log.write("The longest sequence is {0} aa long.\n The shortest sequence is {1} aa long.\n".format(max_seq, min_seq))

### import a list of positions to be used for the alignment; if none given/found: select to align the entire length of
### each sequence or enter specific locations between which sequences should be aligned (slice)
positions = []
try:
    with open(positions_file, "r") as file_object:
        lines = file_object.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            positions.append(line.split(","))
    if positions == []:
        print("No positions could be found.")
        log.write("No positions could be found. User intervention required.\n")
        answer1 = input("Press 'Enter' to align whole sequences or enter a number to indicate tthe start position:")
        if answer1 != "":
            answer2 = input("Enter a number to indicate the end position:")
        else:
            answer2 = ""

except IOError:
    print("'Positions' file not found.")
    log.write("'Positions' file not found.\n")
    answer1 = input("Press 'Enter' to align whole sequences or enter a number to indicate the start position:")
    if answer1 != "":
        answer2 = input("Enter a number to indicate the end position:")
    else:
        answer2 = ""

if answer1 == "":
    for i in range(len(ids)):
        positions.append([1, len(seqs[i])])
    log.write("Entire length of every sequence was chosen to be aligned.")
else:
    for i in range(len(ids)):
        positions.append([int(answer1), int(answer2)])
    log.write("\nSequences will be sliced up between positions {0} and {1}.".format(int(answer1), int(answer2)))

print("List of positions has {0} entries.\n".format(len(positions)))
log.write("\nList of positions has {0} entries.\n".format(len(positions)))

### make sure that the length of positions list equals the number of sequences to be aligned
if len(ids) == len(seqs):
    print("Found and processed {0} entries.".format(len(ids)))
    log.write("Found and processed {0} entries.\n".format(len(ids)))
else:
    print("IDs don't match sequences. {0} IDs and {1} sequences were found.".format(len(ids), len(seqs)))
    log.write("IDs don't match sequences. {0} IDs and {1} sequences were found.\n".format(len(ids), len(seqs)))

### generate a new list containing sequences truncated using the positions given above
trunc_seqs = []
seq_no = len(seqs)
log.write("\n")

for i in range(seq_no):
    log.write("Position {0}:\t{1}\t{2}.\tSequence length: {3}.\tSlice length: {4}.\n".format(i + 1, positions[i][0],
                                                                                             positions[i][1],
                                                                                             len(seqs[i]),
                                                                                             int(positions[i][1]) - int(
                                                                                                 positions[i][0])))

    try:
        pos0 = int(positions[i][0].strip()) - 1
        pos1 = int(positions[i][1].strip())
    except AttributeError:
        pos0 = int(positions[i][0]) - 1
        pos1 = int(positions[i][1])

    if pos1 <= len(seqs[i]):
        trunc_seqs.append(seqs[i][pos0:pos1])
    else:
        trunc_seqs.append(seqs[i])
        print("Sequence {0} is too short for the required truncation. Appending entire sequence.".format(ids[i]))
        log.write(
            "Sequence {0} is too short for the required truncation {1}. Appending entire sequence.\n".format(ids[i],
                                                                                                             positions[
                                                                                                                 i]))

### write a new fasta file containing the sliced sequences - potentially, replace sequence ids with original ids
with open(all_seq_sliced, "w") as file_object:
    for i in range(len(trunc_seqs)):
        file_object.write(">" + str(ids[i]) + "\n" + str(trunc_seqs[i]) + "\n\n")

### use ClustalO to align the sliced sequences
print("\nAligning...")

clustalo_cline = ClustalOmegaCommandline(clustalo_exe, infile=all_seq_sliced, outfile=clustal_out, force=forceO,
                                         verbose=True, auto=True)

assert os.path.isfile(clustalo_exe), "Clustal O executable could not be found" # make sure that ClustalO executable can be found
stdout, stderr = clustalo_cline()

log.write("\nClustalO alignment successfully completed.")

### convert ALN alignment file into FASTA format
with open(clustal_out, 'r') as aln_file:
    aln_fasta_file = open(clustal_out_fasta, 'w')
    aln_fasta_full = open(subfolder + "/full_alignment.fasta", 'w')

    for line in aln_file.readlines():
        aln_fasta_file.write(line)
        aln_fasta_full.write(line)

    aln_fasta_file.close()
    aln_fasta_full.close()

print("Done in {0:.0f} min {1:.0f} s.".format(round((time() - time0) / 60, 0), round(((time() - time0) % 60), 0)))
