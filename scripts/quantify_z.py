from sklearn import svm
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import compartment_analysis as ca
from matplotlib import pyplot as plt
import os
import linear_algebra as la
from scipy import stats as st

res_kb = 100
chroms = (22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 8, 7, 6, 5, 4, 3, 2, 1)
cell_type1 = "GM12878_combined"
cell_type2 = "K562"

x_means = np.zeros_like(chroms, dtype=float)
y_means = np.zeros_like(chroms, dtype=float)
z_means = np.zeros_like(chroms, dtype=float)
x_lengths = np.zeros_like(chroms, dtype=float)
y_lengths = np.zeros_like(chroms, dtype=float)
z_lengths = np.zeros_like(chroms, dtype=float)

for i, chrom in enumerate(chroms):
	path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

	os.system("python ../multimds.py --full {} {}".format(path1, path2))
	structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))
	structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

	#compartments
	contacts1 = dt.matFromBed(path1, structure1)
	contacts2 = dt.matFromBed(path2, structure2)

	compartments1 = np.array(ca.get_compartments(contacts1))
	compartments2 = np.array(ca.get_compartments(contacts2))

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
	coef = clf.coef_

	transformed_coords1 = np.array(la.change_coordinate_system(coef, coords1))
	transformed_coords2 = np.array(la.change_coordinate_system(coef, coords2))

	x_diffs = transformed_coords1[:,0] - transformed_coords2[:,0]
	y_diffs = transformed_coords1[:,1] - transformed_coords2[:,1]
	z_diffs = transformed_coords1[:,2] - transformed_coords2[:,2]

	x_means[i] = np.mean(np.abs(x_diffs))
	y_means[i] = np.mean(np.abs(y_diffs))
	z_means[i] = np.mean(np.abs(z_diffs))

	#axis lengths
	centroid1 = np.mean(transformed_coords1, axis=0)
	centroid2 = np.mean(transformed_coords2, axis=0)
	x_length1 = np.mean([np.abs(coord1[0] - centroid1[0]) for coord1 in transformed_coords1])
	y_length1 = np.mean([np.abs(coord1[1] - centroid1[1]) for coord1 in transformed_coords1])
	z_length1 = np.mean([np.abs(coord1[2] - centroid1[2]) for coord1 in transformed_coords1])
	x_length2 = np.mean([np.abs(coord2[0] - centroid2[0]) for coord2 in transformed_coords2])
	y_length2 = np.mean([np.abs(coord2[1] - centroid2[1]) for coord2 in transformed_coords2])
	z_length2 = np.mean([np.abs(coord2[2] - centroid2[2]) for coord2 in transformed_coords2])

	x_lengths[i] = np.mean((x_length1, x_length2))
	y_lengths[i] = np.mean((y_length1, y_length2))
	z_lengths[i] = np.mean((z_length1, z_length2))

z_fractions = np.zeros_like(z_means)
for i, (x_mean, y_mean, z_mean) in enumerate(zip(x_means, y_means, z_means)):
	z_fractions[i] = z_mean/(x_mean + y_mean + z_mean)

print np.mean(z_fractions)

ind = np.arange(len(chroms))  # the x locations for the groups
width = 0.2       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, x_means, width, color="r")
rects2 = ax.bar(ind + width, y_means, width, color="y")
rects3 = ax.bar(ind + 2*width, z_means, width, color="b")
ax.set_ylabel("Mean absolute change")
ax.set_xticks(ind + width / 3)
ax.set_xticklabels(chroms)
ax.legend((rects1[0], rects2[0], rects3[0]), ("x", "y", "z"))
plt.savefig("change_by_axis")

fig, ax = plt.subplots()
rects1 = ax.bar(ind, x_lengths, width, color="r")
rects2 = ax.bar(ind + width, y_lengths, width, color="y")
rects3 = ax.bar(ind + 2*width, z_lengths, width, color="b")
ax.set_ylabel("Mean length")
ax.set_xticks(ind + width / 3)
ax.set_xticklabels(chroms)
ax.legend((rects1[0], rects2[0], rects3[0]), ("x", "y", "z"))
plt.savefig("axis_length")
