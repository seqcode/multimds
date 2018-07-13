import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
from sklearn.manifold import MDS
import linear_algebra as la
from matplotlib import pyplot as plt
from scipy import stats as st
from joint_mds import Joint_MDS

chrom = 21
res_kb = 100

exp_names = ("GM12878_combined", "K562")

path1 = "/data/drive1/hic_data/{}_{}_{}kb.bed".format(exp_names[0], chrom, res_kb)
path2 = "/data/drive1/hic_data/{}_{}_{}kb.bed".format(exp_names[1], chrom, res_kb)

structure1 = dt.structureFromBed(path1, None, None)
structure2 = dt.structureFromBed(path2, None, None)

dt.make_compatible((structure1, structure2))

#get distance matrices
dists1 = dt.normalized_dist_mat(path1, structure1)
dists2 = dt.normalized_dist_mat(path2, structure2)

n = 10
all_r_sq = []

ps = np.arange(0.01, 0.1, 0.01)
for p in ps:
	all_changes = []
	for i in range(n):
		print i
		#perform MDS
		coords1, coords2 = Joint_MDS(n_components=3, p=p, random_state1=np.random.RandomState(), random_state2=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=-1).fit_transform(dists1, dists2)

		#fill structures
		structure1.setCoords(coords1)
		structure2.setCoords(coords2)

		all_changes.append(np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(structure1.getCoords(), structure2.getCoords())]))

	multimds_r_sq = []
	for i in range(n):
		for j in range(i):
			r, p = st.pearsonr(all_changes[i], all_changes[j])
			multimds_r_sq.append(r**2)

	all_r_sq.append(multimds_r_sq)

all_changes = []
for i in range(n):
	print i
	#perform MDS
	coords1 = MDS(n_components=3, random_state=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=-1).fit_transform(dists1)
	coords2 = MDS(n_components=3, random_state=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=-1).fit_transform(dists2)

	#fill structures
	structure1.setCoords(coords1)
	structure2.setCoords(coords2)

	#kabsch
	r, t = la.getTransformation(structure1, structure2)
	structure1.transform(r,t)

	all_changes.append(np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(structure1.getCoords(), structure2.getCoords())]))

kabsch_r_sq = []
for i in range(n):
	for j in range(i):
		r, p = st.pearsonr(all_changes[i], all_changes[j])
		kabsch_r_sq.append(r**2)

labels = ["p={}".format(p) for p in ps]
labels.append("Kabsch")

all_r_sq.append(kabsch_r_sq)

fig, ax = plt.subplots()
ax.boxplot(all_r_sq)
ax.set_ylabel("Reproducibility", fontsize=20)
ax.set_xticklabels(labels, fontsize=15)
plt.savefig("reproducibility")
