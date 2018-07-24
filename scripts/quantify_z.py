from sklearn import svm
import numpy as np
import data_tools as dt
import sys
sys.path.append("..")
import compartment_analysis as ca
from matplotlib import pyplot as plt
import os
import linear_algebra as la
import array_tools as at
from scipy import stats as st

res_kb = 100
chroms = (22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 8, 7, 6, 5, 4, 3, 2, 1)
cell_type1 = "GM12878_combined"
cell_type2 = "K562"
n = 5

x_corrs = np.zeros_like(chroms, dtype=float)
y_corrs = np.zeros_like(chroms, dtype=float)
z_corrs = np.zeros_like(chroms, dtype=float)
x_means = np.zeros_like(chroms, dtype=float)
y_means = np.zeros_like(chroms, dtype=float)
z_means = np.zeros_like(chroms, dtype=float)
x_lengths = np.zeros_like(chroms, dtype=float)
y_lengths = np.zeros_like(chroms, dtype=float)
z_lengths = np.zeros_like(chroms, dtype=float)

for i, chrom in enumerate(chroms):
	path1 = "/data/drive1/hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "/data/drive1/hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

	min_error = sys.float_info.max
	for iteration in range(n):
		os.system("python ../minimds.py -o {}_ {} {}".format(iteration, path1, path2))

		#load structures
		structure1 = dt.structure_from_file("{}_{}_{}_{}kb_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
		structure2 = dt.structure_from_file("{}_{}_{}_{}kb_structure.tsv".format(iteration, cell_type2, chrom, res_kb))

		#rescale
		structure1.rescale()
		structure2.rescale()

		#make structures compatible
		dt.make_compatible((structure1, structure2))

		#align
		r, t = la.getTransformation(structure1, structure2)
		structure1.transform(r,t)

		#calculate error
		coords1 = np.array(structure1.getCoords())
		coords2 = np.array(structure2.getCoords())
		error = np.mean([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)])
		if error < min_error:
			min_error = error
			best_iteration = iteration

	for iteration in range(n):
		if iteration == best_iteration:
			#load structures
			structure1 = dt.structure_from_file("{}_{}_{}_{}kb_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
			structure2 = dt.structure_from_file("{}_{}_{}_{}kb_structure.tsv".format(iteration, cell_type2, chrom, res_kb))
		else:
			os.system("rm {}_{}_{}_{}kb_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
			os.system("rm {}_{}_{}_{}kb_structure.tsv".format(iteration, cell_type2, chrom, res_kb))	

	#rescale
	structure1.rescale()
	structure2.rescale()

	#make structures compatible
	dt.make_compatible((structure1, structure2))

	#align
	r, t = la.getTransformation(structure1, structure2)
	structure1.transform(r,t)

	#calculate error
	coords1 = np.array(structure1.getCoords())
	coords2 = np.array(structure2.getCoords())
	dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)]

	#compartments
	contacts1 = dt.matFromBed(path1, structure1)
	contacts2 = dt.matFromBed(path2, structure2)

	at.makeSymmetric(contacts1)
	at.makeSymmetric(contacts2)

	compartments1 = ca.infer_compartments(contacts1, structure1, cell_type1, chrom, res_kb)
	compartments2 = ca.infer_compartments(contacts2, structure2, cell_type2, chrom, res_kb)

	#SVR
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

	#correlations
	dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)]

	r, p = st.pearsonr(np.abs(x_diffs), dists)
	x_corrs[i] = r
	r, p = st.pearsonr(np.abs(y_diffs), dists)
	y_corrs[i] = r
	r, p = st.pearsonr(np.abs(z_diffs), dists)
	z_corrs[i] = r

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

fig, ax = plt.subplots()
rects1 = ax.bar(ind, x_corrs, width, color="r")
rects2 = ax.bar(ind + width, y_corrs, width, color="y")
rects3 = ax.bar(ind + 2*width, z_corrs, width, color="b")
ax.set_ylabel("Pearson r")
ax.set_xticks(ind + width / 3)
ax.set_xticklabels(chroms)
ax.legend((rects1[0], rects2[0], rects3[0]), ("x", "y", "z"))
plt.savefig("correlation_by_axis")
plt.show()
