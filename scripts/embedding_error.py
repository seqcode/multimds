import sys
sys.path.append("..")
import data_tools as dt
import array_tools as at
import numpy as np
import linear_algebra as la
from joint_mds import Joint_MDS
from matplotlib import pyplot as plt

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

structure1 = dt.structureFromBed(path1, None, None)
structure2 = dt.structureFromBed(path2, None, None)

#make structures compatible
dt.make_compatible((structure1, structure2))

#get distance matrices
dists1 = dt.normalized_dist_mat(path1, structure1)
dists2 = dt.normalized_dist_mat(path2, structure2)

ps = np.arange(0, 0.6, 0.1)
errors1 = np.zeros_like(ps)
errors2 = np.zeros_like(ps)

for i, p in enumerate(ps):
	#joint MDS
	coords1, coords2 = Joint_MDS(n_components=3, p=p, random_state1=np.random.RandomState(), random_state2=np.random.RandomState(), dissimilarity="precomputed", n_jobs=-1).fit_transform(dists1, dists2)

	errors1[i] = error(dists1, coords1)
	errors2[i] = error(dists2, coords2)

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
