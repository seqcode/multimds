from matplotlib import pyplot as plt
from scipy import stats as st
import numpy as np
import sys

datasets = ("GM12878_H3K27ac", "GM12878_H3K4me1", "GM12878_H3K4me3", "GM12878_H3K9ac", "GM12878_H2AZ", "GM12878_H3K4me2", "Gm12878_enhancer", "K562_H3K27me3", "K562_EZH2", "K562_repressed")

overrepresentation = np.zeros_like(datasets, dtype=float)


for i, dataset in enumerate(datasets):
	peaks_coverage = np.loadtxt("peaks_filtered_independent_{}_coverage.bed".format(dataset), usecols=6)
	background_coverage = np.loadtxt("A_background_independent_filtered_{}_coverage.bed".format(dataset), usecols=6)

	print st.ttest_ind(peaks_coverage, background_coverage)

	overrepresentation[i] = np.mean(peaks_coverage)/np.mean(background_coverage)

ys = overrepresentation 

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 7, 10, frameon=False)

#label axes
plt.ylabel("Relative enrichment", fontsize=12)

#define offsets
xs = range(len(ys))

xmin = min(xs)
xmax = max(xs)
x_range = xmax - xmin
x_start = xmin - x_range/15.	#bigger offset for bar plot
x_end = xmax + x_range/15.

ymin = 0
ymax = max(ys)
y_range = ymax - ymin
y_start = ymin - y_range/50.
y_end = ymax

#plot data
plt.bar(xs, ys, width=0.4, bottom=y_start)

plt.plot([x_start, x_end], [1, 1], linestyle="--", c="k")

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.xticks(xs, ("GM12878 H3K27ac", "GM12878 H3K4me1", "GM12878 H3K4me3", "GM12878 H3K9ac", "GM12878 H2AZ", "GM12878 H3K4me2", "GM12878 enhancer", "K562 H3K27me3", "K562 EZH2", "K562 repressed"), rotation=90)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, labelsize=9)

plt.savefig("sup18")
