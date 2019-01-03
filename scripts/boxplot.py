from matplotlib import pyplot as plt
import sys
import numpy as np

gene_name = sys.argv[1]
sgd = sys.argv[2]

with open("deseq_counts.tsv") as infile:
	for line in infile:
		line = line.strip().split()
		if line[0] == sgd:
			glucose_counts = [np.log10(float(line[i])) for i in range(1,4)]
			galactose_counts = [np.log10(float(line[i])) for i in range(4,7)]

ys = [glucose_counts, galactose_counts]
plt.subplot2grid((10,10), (0,0), 10, 10, frameon=False)

plt.title(gene_name)

#label axes
plt.ylabel("log RNA-seq read count", fontsize=12)

#define offsets
x_start = 0.3
x_end = 0.7

ymin = min([min(y) for y in ys])
ymax = max([max(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, positions=(0.4, 0.6), patch_artist=True, labels=("glucose", "galactose"), medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)

plt.savefig("{}_expression".format(gene_name))
plt.show()	
