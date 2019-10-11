import sys
from multimds import data_tools as dt
from multimds import compartment_analysis as ca
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
import os
import numpy as np

all_species = ("Mouse", "Human", "Yeast")
all_res_kb = (100, 100, 32)

boxes = [[] for species in all_species]

for i, (species, res_kb) in enumerate(zip(all_species, all_res_kb)):
	with open("{}_list.txt".format(species)) as infile:
		for line in infile:
			prefix = line.strip()
			for chrom in range(1, 23):
				path = "hic_data/{}_{}_{}kb.bed".format(prefix, chrom, res_kb)

				if os.path.isfile(path):	
					mat = dt.matFromBed(path)
					oe_mat = ca.oe(mat)
					cor_mat = ca.cor(oe_mat)
					pca = PCA(n_components=1)
					pca.fit(cor_mat)
					boxes[i].append(pca.explained_variance_ratio_[0])

		infile.close()

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("PC1 explained variance ratio", fontsize=10)

#define offsets
ys = boxes
n = len(ys)
width = 0.075

xmin = 0	#boxplot indexing starts at 1
xmax = n*width*2
x_range = xmax - xmin
x_start = xmin - x_range/10.	#larger offset for boxplot
x_end = xmax + x_range/10.

ymin = min([min(y) for y in ys])
ymax = max([max(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, notch=True, patch_artist=True, positions=np.arange(width, n*width*2, width*2), widths=[width for i in range(n)], labels=all_species, medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.savefig("sup9a.svg")
