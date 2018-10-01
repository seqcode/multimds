import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
import linear_algebra as la
from matplotlib import pyplot as plt
from scipy import stats as st
import os

chrom = 21
res_kb = 100
exp_names = ("GM12878_combined", "K562")
path1 = "hic_data/{}_{}_{}kb.bed".format(exp_names[0], chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(exp_names[1], chrom, res_kb)
n = 10

all_r_sq = []

ps = np.arange(0.01, 0.1, 0.01)
for p in ps:
	all_changes = []
	for i in range(n):
		print(i)
		#perform MDS
		os.system("python ../multimds.py --full {} {}".format(path1, path2))

		structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(exp_names[0], chrom, res_kb))
		structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(exp_names[1], chrom, res_kb))

		all_changes.append(np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(structure1.getCoords(), structure2.getCoords())]))

	multimds_r_sq = []
	for i in range(n):
		for j in range(i):
			r, p = st.pearsonr(all_changes[i], all_changes[j])
			multimds_r_sq.append(r**2)

	all_r_sq.append(multimds_r_sq)

all_changes = []
for i in range(n):
	print(i)
	os.system("python minimds.py {}".format(path1))
	os.system("python minimds.py {}".format(path2))

	structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(exp_names[0], chrom, res_kb))
	structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(exp_names[1], chrom, res_kb))

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

fig = plt.figure()  # define the figure window
ax = fig.add_subplot(111)   # define the axis
ax.boxplot(all_r_sq)
ax.set_ylabel("Reproducibility", fontsize=20)
ax.set_xticklabels(labels, fontsize=15)
plt.savefig("reproducibility")
