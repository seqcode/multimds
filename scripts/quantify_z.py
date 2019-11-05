from sklearn import svm
import numpy as np
import sys
from multimds import data_tools as dt
from multimds import compartment_analysis as ca
from matplotlib import pyplot as plt
import os
from multimds import linear_algebra as la
from scipy import stats as st
from multimds import multimds

res_kb = 100
chroms = range(1, int(sys.argv[1]))
design_file = sys.argv[2]
penalty = float(sys.argv[3])

x_means = []
y_means = []
z_means = []
x_lengths = []
y_lengths = []
z_lengths = []

with open(design_file) as infile:
	for line in infile:
		cell_type1, cell_type2 = line.strip().split()
		for chrom in chroms:
			path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
			path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

			if os.path.isfile(path1) and os.path.isfile(path2):
				structure1, structure2 = multimds.full_mds(path1, path2, penalty=penalty)

				structure1.rescale()
				structure2.rescale()
				r,t = la.getTransformation(structure1, structure2)
				structure1.transform(r,t)

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
				coef = clf.coef_

				transformed_coords1 = np.array(la.change_coordinate_system(coef, coords1))
				transformed_coords2 = np.array(la.change_coordinate_system(coef, coords2))

				x_diffs = transformed_coords1[:,0] - transformed_coords2[:,0]
				y_diffs = transformed_coords1[:,1] - transformed_coords2[:,1]
				z_diffs = transformed_coords1[:,2] - transformed_coords2[:,2]

				#axis lengths
				centroid1 = np.mean(transformed_coords1, axis=0)
				centroid2 = np.mean(transformed_coords2, axis=0)
				x_length1 = np.mean([np.abs(coord1[0] - centroid1[0]) for coord1 in transformed_coords1])
				y_length1 = np.mean([np.abs(coord1[1] - centroid1[1]) for coord1 in transformed_coords1])
				z_length1 = np.mean([np.abs(coord1[2] - centroid1[2]) for coord1 in transformed_coords1])
				x_length2 = np.mean([np.abs(coord2[0] - centroid2[0]) for coord2 in transformed_coords2])
				y_length2 = np.mean([np.abs(coord2[1] - centroid2[1]) for coord2 in transformed_coords2])
				z_length2 = np.mean([np.abs(coord2[2] - centroid2[2]) for coord2 in transformed_coords2])

				x_length = np.mean((x_length1, x_length2))
				y_length = np.mean((y_length1, y_length2))
				z_length = np.mean((z_length1, z_length2))

				x_means.append(np.mean(np.abs(x_diffs))/x_length)
				y_means.append(np.mean(np.abs(y_diffs))/y_length)
				z_means.append(np.mean(np.abs(z_diffs))/z_length)


x_fractions = np.zeros_like(x_means)
y_fractions = np.zeros_like(y_means)
z_fractions = np.zeros_like(z_means)
for i, (x_mean, y_mean, z_mean) in enumerate(zip(x_means, y_means, z_means)):
	tot = x_mean + y_mean + z_mean
	x_fractions[i] = x_mean/tot
	y_fractions[i] = y_mean/tot
	z_fractions[i] = z_mean/tot

print st.ttest_ind(x_fractions, y_fractions)
print st.ttest_ind(x_fractions, z_fractions)
print st.ttest_ind(y_fractions, z_fractions)

medianprops = dict(linestyle="none")
labels = ("Orthogonal 1", "Orthogonal 2", "Compartment axis")
prefix = design_file.split("_design.txt")[0]

y_int_size = 0.02
x_start = 0
x_end = 0.6
ymin = 0
ymax = float(sys.argv[4])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.boxplot([x_fractions, y_fractions, z_fractions], positions=(0.075, 0.275, 0.475), widths=(0.1, 0.1, 0.1), notch=True, patch_artist=True, labels=labels, medianprops=medianprops)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.locator_params(axis="y", nbins=2)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)
plt.title(" ".join(sys.argv[5:len(sys.argv)]))
plt.savefig("{}_change_by_axis".format(prefix))
