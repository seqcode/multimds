import sys
sys.path.append("..")
from multimds import linear_algebra as la
from multimds import data_tools as dt
from matplotlib import pyplot as plt
import os
import numpy as np
from ultimds import multimds as mm

def rmsd(struct1, struct2):
	coords1 = struct1.getCoords()
	coords2 = struct2.getCoords()
	assert len(coords1) == len(coords2)
	ssd = 0
	ssd = sum([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)])
	msd = ssd/len(coords1)
	return msd**(1./2)

os.system("cp hic_data/GM12878_combined_21_100kb.bed hic_data/copy-GM12878_combined_21_100kb.bed")

num_iterations = 10
independent_rmsds = np.zeros(num_iterations)
multimds_rmsds = np.zeros(num_iterations)

for i in range(num_iterations):
	os.system("python ../minimds.py hic_data/GM12878_combined_21_100kb.bed")
	os.system("python ../minimds.py hic_data/copy-GM12878_combined_21_100kb.bed")
	struct1 = dt.structure_from_file("hic_data/GM12878_combined_21_100kb_structure.tsv")
	struct2 = dt.structure_from_file("hic_data/copy-GM12878_combined_21_100kb_structure.tsv")
	r,t = la.getTransformation(struct1, struct2)
	struct1.transform(r,t)
	independent_rmsds[i] = rmsd(struct1, struct2)
	
	#os.system("python ../multimds.py -P 0.02 hic_data/GM12878_combined_21_100kb.bed hic_data/copy-GM12878_combined_21_100kb.bed")
	struct1, struct2 = mm.full_mds("hic_data/GM12878_combined_21_100kb.bed", "hic_data/copy-GM12878_combined_21_100kb.bed", penalty=0.02)
	#struct1 = dt.structure_from_file("GM12878_combined_21_100kb_structure.tsv")
	#struct2 = dt.structure_from_file("copy-GM12878_combined_21_100kb_structure.tsv")
	multimds_rmsds[i] = rmsd(struct1, struct2)

ys = [independent_rmsds, multimds_rmsds]

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("RMSD of identical datasets", fontsize=14)

x_start = 0
x_end = 0.4

ymax = max([max(y) for y in ys])
ymin = min([min(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, notch=True, patch_artist=True, labels=("Independent MDS", "MultiMDS"), medianprops=dict(linestyle="none"), positions=(0.075, 0.275), widths=(0.1, 0.1))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)
plt.savefig("sup1")
plt.show()	
