import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
import linear_algebra as la
from matplotlib import pyplot as plt
import os

def error(dists, coords):
	assert len(dists) == len(coords)
	n = len(dists)

	sse = 0
	count = 0
	for i in range(n):
		for j in range(i):
			embedded_dist = la.calcDistance(coords[i], coords[j])
			sse += (embedded_dist - dists[i,j])**2
			count += 1
	mse = sse/count
	rmse = mse**(1./2)
	return rmse

chrom = 21
res_kb = 100
prefix1 = "GM12878_combined"
prefix2 = "K562"
	
path1 = "hic_data/{}_{}_{}kb.bed".format(prefix1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(prefix2, chrom, res_kb)

ps = np.arange(0, 0.6, 0.1)
errors = np.zeros_like(ps)

for i, p in enumerate(ps):
	os.system("python ../multimds.py -P {} {} {}".format(p, path1, path2))
	structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(prefix1, chrom, res_kb))
	structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(prefix2, chrom, res_kb))

	dists1 = dt.normalized_dist_mat(path1, structure1)
	dists2 = dt.normalized_dist_mat(path2, structure2)

	errors[i] = np.mean((error(dists1, structure1.getCoords()), error(dists1, structure1.getCoords())))


xs = ps
x_int_size = 0.1
ys = errors
y_int_size = 0.05
x_start = min(xs) - x_int_size/4.
x_end = max(xs) + x_int_size/5.
y_start = -y_int_size/5.
y_end = max(ys) + y_int_size/5.

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.bar(ps, errors, 0.04, bottom=y_start)
plt.ylabel("Average RMSE", fontsize=14)
plt.xlabel("Similarity weight", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end],frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=-0.01, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)
plt.savefig("{}_{}_embedding_error".format(prefix1, prefix2))
