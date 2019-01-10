import sys
sys.path.append("..")
import compartment_analysis as ca
import data_tools as dt
from scipy import stats as st
import numpy as np
from matplotlib import pyplot as plt
import os

comparisons = ("mouse_celltype", "encode", "cohesin", "lymphoblastoid")
boxes = [[] for comparison in comparisons]

for i, comparison in enumerate(comparisons):
	with open("{}_design.txt".format(comparison)) as infile:
		for line in infile:
			prefix1, prefix2 = line.strip().split()
			for chrom in range(1, 23):
				path1 = "hic_data/{}_{}_100kb.bed".format(prefix1, chrom)
				path2 = "hic_data/{}_{}_100kb.bed".format(prefix2, chrom)

				if os.path.isfile(path1) and os.path.isfile(path2):
					struct1 = dt.structureFromBed(path1)
					struct2 = dt.structureFromBed(path2)

					dt.make_compatible((struct1, struct2))

					mat1 = dt.matFromBed(path1, struct1)
					comps1 = ca.get_compartments(mat1)
					mat2 = dt.matFromBed(path2, struct2)
					comps2 = ca.get_compartments(mat2)

					r, p = st.pearsonr(comps1, comps2)
					if r < 0:
						comps1 = -comps1

					switch = 0
					same = 0
					for comp1, comp2 in zip(comps1, comps2):
						if comp1*comp2 > 0:
							same += np.abs(comp1 - comp2)
						elif comp1*comp2 < 0:
							switch += np.abs(comp1 - comp2)

					boxes[i].append(same/(switch+same))

		infile.close()

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("Intra-compartment change fraction", fontsize=12)

#define offsets
ys = boxes
xs = range(len(ys))

xmin = 1	#boxplot indexing starts at 1
xmax = len(ys)
x_range = xmax - xmin
x_start = xmin - x_range/10.	#larger offset for boxplot
x_end = xmax + x_range/10.

ymin = min([min(y) for y in ys])
ymax = max([max(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, notch=True, patch_artist=True, labels=("Mouse cell types", "ENCODE", "Cohesin KO", "LCLs"), medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.savefig("compartment_strength")
plt.show()
