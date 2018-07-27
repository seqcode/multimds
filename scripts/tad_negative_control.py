import numpy as np
import os
from matplotlib import pyplot as plt
import sys

mat = np.loadtxt("A_background_filtered.bed", dtype=object)
m = len(mat)
ns = []
num_peaks = int(sys.argv[1])
num_overlap = int(sys.argv[2])

for i in range(100):
	indices = np.random.randint(0, m-1, num_peaks)
	rand_mat = mat[indices]
	np.savetxt("negative_control.bed", rand_mat, fmt="%s", delimiter="\t")
	os.system("bedtools intersect -a negative_control.bed -b GM12878_combined_K562_100kb_differential_tad_boundaries.bed > intersection.bed")
	intersection = np.loadtxt("intersection.bed", dtype=object)
	ns.append(len(intersection)/float(num_peaks))

plt.boxplot([ns, [num_overlap/float(num_peaks)]], labels=["Random A compartment", "Relocalization peaks"])
plt.ylabel("Fraction overlap with differential TAD boundaries")
plt.savefig("differential_tad_boundaries_enrichment")
