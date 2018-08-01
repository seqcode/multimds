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
partition_nums = (4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2)
midpoints = (135, 93, 92, 51, 48, 60, 60, 45, 41, 53, 36, 0, 0, 0, 40, 24, 17, 26, 28, 0, 0)
multimds_z_rs = np.zeros_like(chroms, dtype=float)
contacts_pearson_rs = np.zeros_like(chroms, dtype=float)
contacts_spearman_rs = np.zeros_like(chroms, dtype=float)
dists_pearson_rs = np.zeros_like(chroms, dtype=float)
dists_spearman_rs = np.zeros_like(chroms, dtype=float)

for j, chrom in enumerate(chroms):	
	path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

	os.system("python ../multimds.py -m {} -N {} -o hic_data/ {} {}".format(midpoints[j]*10**6, partition_nums[j], path1, path2))

	#load structures
	structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))	
	structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

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
	compartments1 = ca.infer_compartments(mat1, structure1, cell_type1, chrom, res_kb)
	compartments2 = ca.infer_compartments(mat2, structure2, cell_type2, chrom, res_kb)
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

ind = np.arange(len(chroms)) 
width = 0.3       

fig, ax = plt.subplots()
rects1 = ax.bar(ind, contacts_pearson_rs, width, color="r")
rects2 = ax.bar(ind + width, contacts_spearman_rs, width, color="c")
rects3 = ax.bar(ind + 2*width, multimds_z_rs, width, color="orange")

ax.set_xlabel("Chromosome")
ax.set_ylabel("Correlation with compartment changes")
ax.set_xticks(ind + width / 3)
ax.set_xticklabels(chroms)
ax.legend((rects1[0], rects2[0], rects3[0]), ("Pearson", "Spearman", "MultiMDS z"))
plt.savefig("dist_vs_compartment")
