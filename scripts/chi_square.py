from matplotlib import pyplot as plt
from scipy import stats as st
import numpy as np

datasets = ("GM12878_H3K27ac", "GM12878_H3K4me1", "GM12878_H3K4me3", "GM12878_H3K9ac", "GM12878_CTCF", "GM12878_H2AZ", "GM12878_H3K4me2", "GM12878_H3K79me2", "GM12878_H4K20me1", "K562_H3K27me3")

overrepresentation = np.zeros_like(datasets, dtype=float)

for i, dataset in enumerate(datasets):
	peaks_coverage = np.loadtxt("peaks_filtered_{}_coverage.bed".format(dataset), usecols=6)
	background_coverage = np.loadtxt("A_background_filtered_{}_coverage.bed".format(dataset), usecols=6)
	peaks_zero_coverage = len(np.where(peaks_coverage == 0)[0])
	peaks_nonzero_coverage = len(np.where(peaks_coverage > 0)[0])
	background_zero_coverage = len(np.where(background_coverage == 0)[0])
	background_nonzero_coverage = len(np.where(background_coverage > 0)[0])

	chi2, p, dof, expected = st.chi2_contingency([[peaks_zero_coverage, background_zero_coverage], [peaks_nonzero_coverage, background_nonzero_coverage]])
	print p

	overrepresentation[i] = (float(peaks_nonzero_coverage)/peaks_zero_coverage)/(float(background_nonzero_coverage)/background_zero_coverage)

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
plt.xticks(xs, datasets, rotation=90)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, labelsize=9)

plt.savefig("relocalization_enrichment")
plt.show()	
