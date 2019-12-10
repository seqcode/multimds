import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as st
import compartment_analysis as ca
import os

comparisons = ("mouse_celltype", "encode", "cohesin", "lymphoblastoid", "mouse_celltype_rep", "encode_rep")
boxes = [[] for comparison in comparisons]

for i, comparison in enumerate(comparisons):
	print comparison
	with open("{}_design.txt".format(comparison)) as infile:
		for line in infile:
			prefix1, prefix2 = line.strip().split()
			for chrom in range(1, 23):
				path1 = "hic_data/{}_{}_100kb.bed".format(prefix1, chrom)
				path2 = "hic_data/{}_{}_100kb.bed".format(prefix2, chrom)

				if os.path.isfile(path1) and os.path.isfile(path2):

					#load structures
					structure1 = dt.structureFromBed(path1)	
					structure2 = dt.structureFromBed(path2)

					dt.make_compatible((structure1, structure2))
	
					mat1 = dt.matFromBed(path1, structure1)
					mat2 = dt.matFromBed(path2, structure2)

					comps1 = ca.get_compartments(mat1, structure1)
					comps2 = ca.get_compartments(mat2, structure2)

					r, p = st.pearsonr(comps1, comps2)

					boxes[i].append(np.abs(r))

		infile.close()

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("Compartment correlation", fontsize=10)

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
plt.boxplot(ys, notch=True, patch_artist=True, positions=np.arange(width, n*(width*2), width*2), widths=[width for i in range(n)], labels=("Mouse cell types", "ENCODE", "Cohesin KO", "LCLs", "Mouse cell type reps", "GM12878 reps"), medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.savefig("sup13")
plt.show()
