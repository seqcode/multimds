import numpy as np
from .linear_algebra import *
from matplotlib import pyplot as plt
from scipy import stats as st
from .multimds import full_mds
import os

def plot_reproducibility(path1, path2):
	prefix1 = os.path.basename(path1.split(".")[0])
	prefix2 = os.path.basename(path2.split(".")[0])

	n = 10

	all_r_sq = []

	ps = np.arange(0, 0.1, 0.01)

	for p in ps:
		all_changes = []
		for i in range(n):
			structure1, structure2 = full_mds(path1, path2, penalty=p)
			
			if p == 0:
				r, t = getTransformation(structure1, structure2)
				structure1.transform(r,t)

			all_changes.append(np.array([calcDistance(coord1, coord2) for coord1, coord2 in zip(structure1.getCoords(), structure2.getCoords())]))

		r_sq = []
		for i in range(n):
			for j in range(i):
				r, p = st.pearsonr(all_changes[i], all_changes[j])
				r_sq.append(r**2)

		all_r_sq.append(r_sq)

	ys = all_r_sq

	#start with a frameless plot (extra room on the left)
	plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

	#label axes
	plt.xlabel("Similarity weight", fontsize=14)
	plt.ylabel("Correlation between iterations", fontsize=14)

	#define offsets
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
	plt.boxplot(ys, notch=True, patch_artist=True, labels=ps, medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

	#define axes with offsets
	plt.axis([x_start, x_end, y_start, y_end], frameon=False)

	#plot axes (black with line width of 4)
	plt.axvline(x=x_start, color="k", lw=4)
	plt.axhline(y=y_start, color="k", lw=4)

	#plot ticks
	plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

	plt.savefig("{}_{}_reproducibility".format(prefix1, prefix2))
