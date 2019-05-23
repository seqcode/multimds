import os
import sys
sys.path.append("..")
import data_tools as dt
import linear_algebra as la
from sklearn import svm
import numpy as np
import compartment_analysis as ca
from scipy import stats as st
from matplotlib import pyplot as plt

res_kb = 100
cell_type1 = "GM12878_combined"
cell_type2 = "K562"
penalty = 0.025

chroms = (21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 8, 7, 6, 5, 4, 3, 2, 1)
multimds_coeffs = np.zeros_like(chroms, dtype=float)
unaligned_coeffs = np.zeros_like(multimds_coeffs)

for i, chrom in enumerate(chroms):

	path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

	os.system("python ../multimds.py -P {} {} {}".format(penalty, path1, path2))
	structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))
	structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))
                             
	#compartments
	contacts1 = dt.matFromBed(path1, structure1)
	contacts2 = dt.matFromBed(path2, structure2)

	compartments1 = np.array(ca.get_compartments(contacts1, 1))
	compartments2 = np.array(ca.get_compartments(contacts2, 1))

	r, p = st.pearsonr(compartments1, compartments2)
	if r < 0:
		compartments2 = -compartments2

	#SVR
	coords1 = structure1.getCoords()
	coords2 = structure2.getCoords()
	coords = np.concatenate((coords1, coords2))
	compartments = np.concatenate((compartments1, compartments2))
	clf = svm.LinearSVR()
	clf.fit(coords, compartments)
	multimds_coeffs[i] = clf.score(coords, compartments)

	os.system("python ../minimds.py {}".format(path1))
	structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))
	os.system("python ../minimds.py {}".format(path2))
	structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

	dt.make_compatible((structure1, structure2))

	#SVR
	coords1 = structure1.getCoords()
	coords2 = structure2.getCoords()
	coords = np.concatenate((coords1, coords2))
	clf = svm.LinearSVR()
	clf.fit(coords, compartments)
	unaligned_coeffs[i] = clf.score(coords, compartments)

print np.mean(multimds_coeffs)
print np.mean(unaligned_coeffs)
ys = [multimds_coeffs, unaligned_coeffs]

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("SVR R^2", fontsize=14)

x_start = 0
x_end = 0.4

ymax = 1
ymin = min([min(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, notch=True, patch_artist=True, labels=("Aligned", "Unaligned"), medianprops=dict(linestyle="none"), positions=(0.075, 0.275), widths=(0.1, 0.1))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.savefig("sup5b")
