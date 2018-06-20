import pandas as pd
import numpy as np
from Bio import SeqIO
from sklearn.decomposition import PCA
from time import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pca_utils import eucl_dist, positions2color

### indicate main subfolders where all the data is/will be stored:
subfolder = "trx_sequences/analysis"
subfolder2 = "trx_sequences"
plots_subfolder = subfolder + "/plots"


time0 = time()

### read in the data
aln_df = pd.read_csv(subfolder + "/ranking_table.csv", sep=',', header=0)
org_data = pd.read_csv(subfolder + "/original_data.csv", sep=',', header=0)

dat = pd.DataFrame(aln_df) # makes a copy of ranking table (?)
n = len(dat.columns) # get the number of positions in the alignment

### perform PCA
pca = PCA(n_components=int(0.1*n)) # only keep the amount of principal components equal to 10% of all positions
pca.fit(dat)

print("PCA finished in {0:.0f} min {1:.0f} s.".format(round((time()-time0)/60,0),round(((time()-time0)%60),0)))
### project data into PC space

# 0,1,3 denote PC1, PC2 and PC3; change values for other PCs
xvector = pca.components_[0]
yvector = pca.components_[1]
zvector = pca.components_[2]
bvector = pca.components_[3]

xs = pca.transform(dat)[:,0]
ys = pca.transform(dat)[:,1]
zs = pca.transform(dat)[:,2]
bs = pca.transform(dat)[:,3]

#################------------- GENERATE FIGURES -------------#################

### FIGURE 1 - PC scores and explained variance

scores_plot_2D = plt.figure(figsize=(15, 13))

scores_ax1 = scores_plot_2D.add_subplot(2, 2, 1) # PC1 vs PC2
scores_ax2 = scores_plot_2D.add_subplot(2, 2, 2) # PC1 vs PC3
scores_ax3 = scores_plot_2D.add_subplot(2, 2, 3) # PC2 vs PC3

cumul_ax = scores_plot_2D.add_subplot(2, 2, 4) # explained variance - cumulative
variance_ax = scores_plot_2D.add_subplot(2, 2, 4) # explained variance

# plot explained variance
expl_var = pca.explained_variance_ratio_
pc2plot = 10 # how many PCs should be plotted

xval_var = range(1, pc2plot + 1) # values for the x-axis
yval_var = expl_var[:pc2plot] # values for the y-axis
yval_cum = np.cumsum(yval_var) # values for the y-axis - cumulative

cumul_ax.plot(xval_var, yval_cum, color='r', label="cumulative")
variance_ax.bar(xval_var, yval_var)

# label axes and legend
cumul_ax.set_xlabel("Component", fontsize=16, fontname="Arial")
cumul_ax.set_ylabel("Explained variance", fontsize=16, fontname="Arial")
cumul_ax.legend(loc="best", fontsize=14, prop={'family': 'cursive', 'weight': 'roman', 'size': 'large'})
cumul_ax.set_xticks(np.arange(1,pc2plot+1))


# plot PCA scores
scores_ax1.scatter(xs, ys, c='r', s=30)
scores_ax2.scatter(xs, zs, c='b', s=30)
scores_ax3.scatter(ys, zs, c='g', s=30)

scores_ax1.set_xlabel("PC1", fontsize=16, fontname="Arial")
scores_ax1.set_ylabel("PC2", fontsize=16, fontname="Arial")
scores_ax2.set_xlabel("PC1", fontsize=16, fontname="Arial")
scores_ax2.set_ylabel("PC3", fontsize=16, fontname="Arial")
scores_ax3.set_xlabel("PC2", fontsize=16, fontname="Arial")
scores_ax3.set_ylabel("PC3", fontsize=16, fontname="Arial")

for ax in (scores_ax1, scores_ax2, scores_ax3, variance_ax):
    for tick_x in ax.xaxis.get_major_ticks():
        tick_x.label.set_fontsize(14)
    for tick_xl in ax.xaxis.get_ticklabels():
        tick_xl.set_fontname('Arial')
    for tick_y in ax.yaxis.get_major_ticks():
        tick_y.label.set_fontsize(14)
    for tick_yl in ax.yaxis.get_ticklabels():
        tick_yl.set_fontname('Arial')

plt.suptitle('PCA scores', fontsize=20)
plt.tight_layout()

plt.show()

# save the figures
scores_plot_2D.savefig(subfolder + "/Fig1_scores.png", dpi=1200)
#scores_plot_2D.savefig(subfolder + "/Fig1_scores.eps", dpi=1200)

### FIGURE 2 - loadings plot in 2D

# calculate loadings
loadings = -1 * pca.components_.T * np.sqrt(pca.explained_variance_)

loadX = []
loadY = []

for i in range(len(loadings)):
    loadX.append(loadings[i][1])  # -----> PCs are swapped here!
    loadY.append(loadings[i][0])

loadings_figure_2D = plt.figure(figsize=(10, 10))

ax1 = loadings_figure_2D.add_subplot(1, 1, 1)
ax1.axhline(0, linestyle='--', color='grey', linewidth=0.5)
ax1.axvline(0, linestyle='--', color='grey', linewidth=0.5)

variable_positions_pc1 = []
nonvariable_positions_pc1 = []
all_positions_X = []
all_positions_X_vals = []
positions_cut1 = []
positions_cons1 = []

# find max euclidean distance in the PC1 space
e_distances = [eucl_dist(x, y) for (x, y) in zip(loadX, loadY)]
e_distances_X = [eucl_dist(x, 0) for x in loadX]

max_e_X = max(e_distances_X) # find the maximum distance to the origin

for i in range(len(loadX)):
    # circles project sequences (ie rows from csv) as points onto PC axes
    x_coord = loadX[i]
    y_coord = loadY[i]

    cutoff_cons1 = 0.02 * max_e_X
    cutoff_var1 = 0.6 * max_e_X

    all_positions_X.append(i)
    all_positions_X_vals.append(x_coord)

    # for eucl distances using only PC1:
    x_coord_abs = np.absolute(x_coord)
    if (x_coord_abs >= cutoff_var1):
        color_load = 'darkmagenta'
        positions_cut1.append([i, x_coord])
    elif (x_coord_abs <= cutoff_cons1):
        color_load = 'dodgerblue'
        positions_cons1.append([i, x_coord])
    else:
        color_load = 'k'

    if (x_coord < -cutoff_var1) or (x_coord > cutoff_var1):
        variable_positions_pc1.append([i, x_coord])
    else:
        nonvariable_positions_pc1.append(i)

    # plot
    ax1.scatter(y_coord, x_coord, c=color_load, s=20)

def relabel_plot2(axis):
    patch1 = mpatches.Patch(color='dodgerblue', label='<2%')
    patch2 = mpatches.Patch(color='darkmagenta', label='>60%')
    axis.legend(loc='upper left', ncol=1, fontsize=16,
                handles=[patch1, patch2],
                prop={'family': 'cursive', 'weight': 'roman', 'size': 'x-large'})


relabel_plot2(ax1)

for tick_x in ax1.xaxis.get_major_ticks():
    tick_x.label.set_fontsize(16)
for tick_xl in ax1.xaxis.get_ticklabels():
    tick_xl.set_fontname('Arial')
for tick_y in ax1.yaxis.get_major_ticks():
    tick_y.label.set_fontsize(16)
for tick_yl in ax1.yaxis.get_ticklabels():
    tick_yl.set_fontname('Arial')

ax1.set_xlabel("PC1", fontsize=18)
ax1.set_ylabel("PC2", fontsize=18)
# ax2.set_xlabel("PC1")
# ax2.set_ylabel("PC3")
# ax3.set_xlabel("PC2")
# ax3.set_ylabel("PC3")

plt.suptitle('Factor loadings', fontsize=20)
plt.show()

loadings_figure_2D.savefig(subfolder + "/Fig2_loadings.png", dpi=1200)
#loadings_figure_2D.savefig(subfolder + "/Fig2_loadings.eps", dpi=1200)

############################### transform the 'alignment positions' into real positions and plot again

def get_positions(list_with_positions, list_with_values):
    l = []
    for i in range(len(list_with_positions)):
        l.append([list_with_positions[i], list_with_values[list_with_positions[i]]])
    return l

def locate_positions(x_coords, y_coords, conserved_coff=0.1, variable_coff=0.4):
    def eucl_dist(xval, yval):
        from numpy import sqrt
        return sqrt(xval ** 2 + yval ** 2)

    e_distances = [eucl_dist(x, y) for (x, y) in zip(x_coords, y_coords)]
    max_e = max(e_distances)
    conserved_positions = []
    variable_positions = []
    other_positions = []

    for i in range(len(x_coords)):
        x = x_coords[i]
        y = y_coords[i]
        e_dist = eucl_dist(x, y)
        if e_dist < conserved_coff * max_e:
            conserved_positions.append([i, x, y])
        elif e_dist > variable_coff * max_e:
            variable_positions.append([i, x, y])
        else:
            other_positions.append([i, x, y])
    return conserved_positions, variable_positions, other_positions


cons_pos_both, var_pos_both, left_pos_both = locate_positions(loadX, loadY, conserved_coff=0.1, variable_coff=0.5)

### convert positions from alignment to absolute positions in a given sequence
seq_id = 197 #which sequence from the data table should be used

record_iterator = SeqIO.parse(subfolder2 + "/trx_sequences_all.fasta", "fasta")
ids = []
seqs = []
for record in record_iterator:
    ids.append(record.id)
    seqs.append(record.seq)

print("Selected sequence: {0}\n{1}".format(ids[seq_id], seqs[seq_id]))

seq_aln = org_data.iloc[seq_id]
seq_raw = seq_aln[seq_aln != '-'] #removes gaps but preserves original indexing
print("Len of aligned seq: {0}\nLen of degapped (original) seq: {1}\n".format(len(seq_aln),len(seq_raw)))

positions = [int(x[3:]) for x in list(seq_raw.index.values)] #contains aa positions from the original sequence (no gaps) numbered using indices from the alignment
gaps = len(seq_aln)-len(seq_raw)

print("No of variable positions in both directions: {0}\nNo of variable positions in PC1 direction: {1}".format(
    len(variable_positions), len(variable_positions_pc1)))

cons_raw = positions2color(cons_pos_both, seq_aln)
var_raw = positions2color(var_pos_both, seq_aln)
rest_raw = positions2color(left_pos_both, seq_aln)


def print_positions(df, label_string):
    print(label_string + " = " + str(list(df['PosRaw'])))


print_positions(cons_raw, "my_list1")
print_positions(var_raw, "my_list2")
