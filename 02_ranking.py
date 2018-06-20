from Bio import SeqIO
import pandas as pd
from time import time
from pca_utils import rank_alphabet, my_split

# indicate main subfolders where all the data is/will be stored:
subfolder = "trx_sequences/analysis"
subfolder2 = "trx_sequences"

time0 = time()

master_list = []
for seq_record in SeqIO.parse(subfolder + "/full_alignment.fasta", "fasta"):
    al_seq = my_split(str(seq_record.seq))
    master_list.append(al_seq)

col_names = [('pos' + str(x)) for x in range(len(master_list[0]))]

aln_df = pd.DataFrame(master_list)  # turn the list of lists into a DataFrame
aln_df.columns = col_names
org_data = aln_df.copy()  # stores a copy of the original dataset

occur_freq = aln_df.apply(pd.value_counts).fillna(0)

################################# replace amino acid letters with corresponding ranks ##############################
print("Ranking positions...")
for column in col_names:
    new_ranks = rank_alphabet(occur_freq[column])
    aln_df[column].replace(new_ranks, inplace=True)

### aln_df is now a DataFrame containing original amino acid sequences replaced with corresponding amino acid ranks

# correct the ranks for the presence of gaps
print("Substituting gaps...")
for j in range(len(col_names)):
    aln_df.loc[org_data[col_names[j]] == '-', col_names[j]] = 0

print("Ranking matrix generated in {0:.0f} min {1:.0f} s.".format(round((time() - time0) / 60, 0),
                                                                  round(((time() - time0) % 60), 0)))

### export the ranking table to csv file

aln_df.to_csv(subfolder + "/ranking_table.csv", sep=',', index=False, header=True)
org_data.to_csv(subfolder + "/original_data.csv", sep=',', index=False, header=True)

print("Table export succesful.")
