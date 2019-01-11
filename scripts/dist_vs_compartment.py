import sys
sys.path.append("..")
from matplotlib import pyplot as plt	
import data_tools as dt
import numpy as np
import compartment_analysis as ca
from scipy import stats as st
import linear_algebra as la
import os
from sklearn import svm

res_kb = 100
cell_type1 = "GM12878_combined"
cell_type2 = "K562"
chroms = (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
multimds_z_rs = np.zeros_like(chroms, dtype=float)
contacts_pearson_rs = np.zeros_like(chroms, dtype=float)
contacts_spearman_rs = np.zeros_like(chroms, dtype=float)

for j, chrom in enumerate(chroms):	
	path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

	#os.system("python ../multimds.py {} {}".format(path1, path2))

	#load structures
	structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))	
	structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

	#compartments
	mat1 = dt.matFromBed(path1, structure1)
	mat2 = dt.matFromBed(path2, structure2)

	compartments1 = ca.get_compartments(mat1)
	compartments2 = ca.get_compartments(mat2)
	r, p = st.pearsonr(compartments1, compartments2)
	if r < 0:
		compartments2 = -compartments2
	compartment_diffs = compartments1 - compartments2

	#SVR
	coords1 = structure1.getCoords()
	coords2 = structure2.getCoords()

	coords = np.concatenate((coords1, coords2))
	compartments = np.concatenate((compartments1, compartments2))
	clf = svm.LinearSVR()
	clf.fit(coords, compartments)
	coef = clf.coef_

	transformed_coords1 = np.array(la.change_coordinate_system(coef, coords1))
	transformed_coords2 = np.array(la.change_coordinate_system(coef, coords2))

	z_diffs = transformed_coords1[:,2] - transformed_coords2[:,2]
	r, p = st.pearsonr(z_diffs, compartment_diffs)
	multimds_z_rs[j] = r

	#contacts Pearson
	rs = np.zeros(len(mat1))
	for i, (row1, row2) in enumerate(zip(mat1, mat2)):
		rs[i], p = st.pearsonr(row1, row2)

	r, p = st.pearsonr(1-rs, np.abs(compartment_diffs))
	contacts_pearson_rs[j] = r

	#contacts Spearman
	rs = np.zeros(len(mat1))
	for i, (row1, row2) in enumerate(zip(mat1, mat2)):
		rs[i], p = st.spearmanr(row1, row2)

	r, p = st.pearsonr(1-rs, np.abs(compartment_diffs))
	contacts_spearman_rs[j] = r

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 10, 10, frameon=False)

#label axes
plt.xlabel("Chromosome", fontsize=13)
plt.ylabel("Correlation with compartment changes", fontsize=13)

#define offsets
xs = np.arange(len(chroms))

xmin = min(xs)
xmax = max(xs)
x_range = xmax - xmin
x_start = xmin - x_range/15.	#bigger offset for bar plot
x_end = xmax + x_range/15.

ymin = min([min(multimds_z_rs), min(contacts_pearson_rs), min(contacts_spearman_rs)])
#ymax = max([max(multimds_z_rs), max(contacts_pearson_rs), max(contacts_spearman_rs)])
ymax = 1
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

width = 0.2

#plot data
plt.bar(xs, multimds_z_rs, width=width, label="MultiMDS")
plt.bar(xs+width, contacts_pearson_rs, width=width, label="Vector pearson r")
plt.bar(xs+2*width, contacts_spearman_rs, width=width, label="Vector spearman r")

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=0, color="k", lw=2)

#plot ticks
plt.xticks(xs, chroms)
plt.tick_params(direction="out", top=False, right=False, bottom=False, length=12, width=3, pad=1, labelsize=10)

plt.legend(frameon=False, fontsize=8, loc=0)

plt.savefig("dist_vs_compartment")
plt.show()
