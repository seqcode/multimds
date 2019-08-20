from matplotlib import pyplot as plt
from scipy import stats as st
import numpy as np

marks = ("H3K27ac", "H3K4me1", "H3K4me3", "H3K9ac", "H3K36me3", "H2AZ", "H3K4me2", "H3K79me2", "H3K27me3", "EZH2", "enhancer", "transcription", "polycomb")
datasets = []
for celltype in ("GM12878", "K562"):
	for mark in marks:
		datasets.append("{}_{}".format(celltype, mark))

overrepresentation = np.zeros_like(datasets, dtype=float)
ps = np.zeros_like(datasets, dtype=float)

for i, dataset in enumerate(datasets):
	peaks_coverage = np.loadtxt("peaks_b_comp_filtered_independent_{}_coverage.bed".format(dataset), usecols=8)
	background_coverage = np.loadtxt("B_background_independent_filtered_{}_coverage.bed".format(dataset), usecols=8)
	
	t, p = st.ttest_ind(peaks_coverage, background_coverage)
	ps[i] = p

	overrepresentation[i] = (np.mean(peaks_coverage)/np.mean(background_coverage) - 1)*100

ys = overrepresentation 

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 7, 10, frameon=False)

#label axes
plt.ylabel("Percent enrichment", fontsize=12)

#define offsets
xs = range(len(ys))

xmin = min(xs)
xmax = max(xs)
x_range = xmax - xmin
x_start = xmin - x_range/15.	#bigger offset for bar plot
x_end = xmax + x_range/15.

ymin = min(ys)
ymax = max(ys)
y_range = ymax - ymin
y_start = ymin - y_range/50.
y_end = ymax + y_range/30.

#plot data
plt.bar(xs, ys, width=0.4, bottom=0)

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=0, color="k", lw=3)

for x, p in zip(xs, ps):
	if p < 0.01:
		plt.scatter([x], [ymax + 0.5], marker="*", s=6, color="k")

#plot ticks
plt.xticks(xs, datasets, rotation=90)
plt.tick_params(direction="out", top=False, right=False, bottom=False, length=12, width=3, labelsize=9)

plt.savefig("sup18_b_comp")
plt.show()
