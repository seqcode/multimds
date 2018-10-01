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
errors1 = np.zeros_like(ps)
errors2 = np.zeros_like(ps)

for i, p in enumerate(ps):
	os.system("python ../multimds.py -P {} --full {} {}".format(p, path1, path2))
	structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(prefix1, chrom, res_kb))
	structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(prefix2, chrom, res_kb))

	dists1 = dt.normalized_dist_mat(path1, structure1)
	dists2 = dt.normalized_dist_mat(path2, structure2)

	errors1[i] = error(dists1, structure1.getCoords())
	errors2[i] = error(dists2, structure2.getCoords())

errors1 = errors1/errors1[0]
errors2 = errors2/errors2[0]

width = 0.35
ind = np.arange(len(ps)) 

fig, ax = plt.subplots()
rects1 = ax.bar(ind, errors1, width, color="b")
rects2 = ax.bar(ind+width, errors2, width, color="g")
plt.ylabel("Normalized RMSE")
plt.ylabel("p")
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(ps)
ax.legend((rects1[0], rects2[0]), ("Embedding 1", "Embedding 2"))
plt.savefig("embedding_error")
