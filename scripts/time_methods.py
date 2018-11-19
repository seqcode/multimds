from matplotlib import pyplot as plt
import numpy as np

chrom_sizes = np.loadtxt("chrom_sizes.txt")

for method in ("multimds", "kabsch"):
	times = np.loadtxt("{}_times.txt".format(method))/60.
	plt.scatter(chrom_sizes, times, label=method)

plt.xlabel("Number of bins")
plt.ylabel("Computational time (minutes)")
plt.legend()
plt.savefig("time_methods")
