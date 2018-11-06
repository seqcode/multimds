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

ps = np.arange(0, 0.1, 0.01)

for p in ps:
	all_changes = []
	for i in range(n):
		print(i)
		os.system("python ../multimds.py -P {} --full {} {}".format(p, path1, path2))

		structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(exp_names[0], chrom, res_kb))
		structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(exp_names[1], chrom, res_kb))
		
		if p == 0:
			r, t = la.getTransformation(structure1, structure2)
			structure1.transform(r,t)

		all_changes.append(np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(structure1.getCoords(), structure2.getCoords())]))

	r_sq = []
	for i in range(n):
		for j in range(i):
			r, p = st.pearsonr(all_changes[i], all_changes[j])
			r_sq.append(r**2)

	all_r_sq.append(r_sq)

medianprops = dict(linestyle="none")
plt.boxplot(all_r_sq, notch=True, patch_artist=True, labels=ps, medianprops=medianprops, showfliers=False)
plt.ylabel("Correlation between iterations", fontsize=15)
plt.xlabel("Difference penalty", fontsize=15)
plt.savefig("reproducibility")
