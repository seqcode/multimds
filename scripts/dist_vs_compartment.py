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
chroms = (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
multimds_z_rs = np.zeros_like(chroms, dtype=float)
contacts_pearson_rs = np.zeros_like(chroms, dtype=float)
contacts_spearman_rs = np.zeros_like(chroms, dtype=float)

for j, chrom in enumerate(chroms):	
	path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

	os.system("python ../multimds.py --full {} {}".format(path1, path2))

	#load structures
	structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))	
	structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

	#rescale
	structure1.rescale()
	structure2.rescale()

	#make structures compatible
	dt.make_compatible((structure1, structure2))

	coords1 = structure1.getCoords()
	coords2 = structure2.getCoords()

	mat1 = dt.matFromBed(path1, structure1)
	mat2 = dt.matFromBed(path2, structure2)
	
	#compartments
	compartments1 = ca.get_compartments(mat1, 1)
	compartments2 = ca.get_compartments(mat2, 1)
	r, p = st.pearsonr(compartments1, compartments2)
	if r < 0:
		compartments2 = -compartments2
	compartment_diffs = compartments1 - compartments2

	#SVR
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

	contacts1 = dt.matFromBed(path1, structure1)
	contacts2 = dt.matFromBed(path2, structure2)

	#contacts Pearson
	rs = np.zeros(len(contacts1))
	for i, (row1, row2) in enumerate(zip(contacts1, contacts2)):
		rs[i], p = st.pearsonr(row1, row2)

	r, p = st.pearsonr(1-rs, np.abs(compartment_diffs))
	contacts_pearson_rs[j] = r

	#contacts Spearman
	rs = np.zeros(len(contacts1))
	for i, (row1, row2) in enumerate(zip(contacts1, contacts2)):
		rs[i], p = st.spearmanr(row1, row2)

	r, p = st.pearsonr(1-rs, np.abs(compartment_diffs))
	contacts_spearman_rs[j] = r

xs = np.arange(len(chroms))
x_int_size = 1
y_int_size = 0.1
x_start = min(xs) - x_int_size/4.
x_end = max(xs) + x_int_size
y_start = -y_int_size/5.
y_end = max((max(multimds_z_rs), max(contacts_pearson_rs), max(contacts_spearman_rs))) + y_int_size/5.

width = 0.2
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
rects1 = plt.bar(xs, multimds_z_rs, width, bottom=y_start, color="g")
rects2 = plt.bar(xs + width, contacts_pearson_rs, width, bottom=y_start, color="b")
rects3 = plt.bar(xs + 2*width, contacts_spearman_rs, width, bottom=y_start, color="orange")
plt.ylabel("Correlation with compartment changes", fontsize=14)
plt.xlabel("Chromosome", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end],frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.xticks(xs, chroms)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)
plt.legend((rects1[0], rects2[0], rects3[0]), ("MultiMDS", "Pearson", "Spearman"))
plt.savefig("dist_vs_compartment")
