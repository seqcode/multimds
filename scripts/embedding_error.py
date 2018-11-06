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

chrom = sys.argv[1]
res_kb = 100
prefix1 = "GM12878_combined"
prefix2 = "K562"
	
path1 = "hic_data/{}_{}_{}kb.bed".format(prefix1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(prefix2, chrom, res_kb)

ps = np.arange(0, 0.6, 0.1)
errors = np.zeros_like(ps)

for i, p in enumerate(ps):
	os.system("python ../multimds.py -P {} --full {} {}".format(p, path1, path2))
	structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(prefix1, chrom, res_kb))
	structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(prefix2, chrom, res_kb))

	dists1 = dt.normalized_dist_mat(path1, structure1)
	dists2 = dt.normalized_dist_mat(path2, structure2)

	errors[i] = np.mean((error(dists1, structure1.getCoords()), error(dists1, structure1.getCoords())))

ind = np.arange(len(ps)) 
width = 0.35

plt.subplot2grid((10,10), (0,0), 10, 10, frameon=False)
plt.bar(ind+width, errors, width)
plt.ylabel("Average RMSE", fontsize=15)
plt.xlabel("Difference penalty", fontsize=15)
plt.axvline(color="k", lw=4)
plt.axhline(color="k", lw=6)
plt.xticks(ind, ps)
plt.savefig("embedding_error")
