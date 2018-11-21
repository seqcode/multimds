import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
import linear_algebra as la
from matplotlib import pyplot as plt
from scipy import stats as st
import os

path1 = sys.argv[1]
path2 = sys.argv[2]

prefix1 = os.path.basename(path1.split(".")[0])
prefix2 = os.path.basename(path2.split(".")[0])

n = 10

all_r_sq = []

ps = np.arange(0, 0.1, 0.01)

for p in ps:
	all_changes = []
	for i in range(n):
		print(i)
		os.system("python ../multimds.py -P {} --full {} {}".format(p, path1, path2))

		structure1 = dt.structure_from_file("{}_structure.tsv".format(os.path.basename(prefix1)))
		structure2 = dt.structure_from_file("{}_structure.tsv".format(os.path.basename(prefix2)))
		
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

x_int_size = 0.01
y_int_size = 0.2
x_start = min(ps) - x_int_size/5.
x_end = max(ps)
y_start = -y_int_size/5.
y_end = 1

medianprops = dict(linestyle="none")
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.boxplot(all_r_sq, notch=True, patch_artist=True, labels=ps, medianprops=medianprops)
plt.ylabel("Correlation between iterations", fontsize=14)
plt.xlabel("Difference penalty", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)
plt.savefig("{}_{}_reproducibility".format(prefix1, prefix2))
plt.show()
