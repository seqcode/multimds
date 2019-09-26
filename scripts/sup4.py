import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

chrom_sizes = np.loadtxt("chrom_sizes.txt")

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
maxs = []
for method in ("MultiMDS", "Independent_MDS"):
	times = np.loadtxt("{}_times.txt".format(method))/60.
	plt.scatter(chrom_sizes, times, label=" ".join(method.split("_")))
	maxs.append(max(times))

xs = chrom_sizes
x_int_size = 500
y_int_size = 1
x_start = min(xs) - x_int_size/4.
x_end = max(xs) + x_int_size/5.
y_start = -y_int_size/5.
y_end = max(maxs) + y_int_size/5.

plt.xlabel("Number of bins", fontsize=14)
plt.ylabel("Computational time (minutes)", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)
plt.legend()
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)
plt.savefig("sup4.svg")
