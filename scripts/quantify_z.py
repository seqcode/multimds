from sklearn import svm
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import compartment_analysis as ca
from matplotlib import pyplot as plt
import os
import linear_algebra as la
import array_tools as at
from scipy import stats as st

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
				os.system("python ../multimds.py --full -P {} {} {}".format(penalty, path1, path2))
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

				#x_lengths.append(np.mean((x_length1, x_length2)))
				#y_lengths.append(np.mean((y_length1, y_length2)))
				#z_lengths.append(np.mean((z_length1, z_length2)))

x_fractions = np.zeros_like(x_means)
y_fractions = np.zeros_like(y_means)
z_fractions = np.zeros_like(z_means)
for i, (x_mean, y_mean, z_mean) in enumerate(zip(x_means, y_means, z_means)):
	tot = x_mean + y_mean + z_mean
	x_fractions[i] = x_mean/tot
	y_fractions[i] = y_mean/tot
	z_fractions[i] = z_mean/tot

print np.mean(z_fractions)
print st.ttest_ind(x_fractions, y_fractions)
print st.ttest_ind(x_fractions, z_fractions)
print st.ttest_ind(y_fractions, z_fractions)

#x_length_fractions = np.zeros_like(x_lengths)
#y_length_fractions = np.zeros_like(y_lengths)
#z_length_fractions = np.zeros_like(z_lengths)
#for i, (x_length, y_length, z_length) in enumerate(zip(x_lengths, y_lengths, z_lengths)):
#	tot = x_length + y_length + z_length
#	x_length_fractions[i] = x_length/tot
#	y_length_fractions[i] = y_length/tot
#	z_length_fractions[i] = z_length/tot

#print np.mean(z_length_fractions)
#print st.ttest_ind(x_length_fractions, y_length_fractions)
#print st.ttest_ind(x_length_fractions, z_length_fractions)
#print st.ttest_ind(y_length_fractions, z_length_fractions)

medianprops = dict(linestyle="none")
labels = ("Orthogonal 1", "Orthogonal 2", "Compartment axis")
prefix = design_file.split("_design.txt")[0]

y_int_size = 0.02
x_start = 0.5
x_end = 3.5
y_start = min((min(x_fractions), min(y_fractions), min(z_fractions))) -y_int_size/5.
y_end = max((max(x_fractions), max(y_fractions), max(z_fractions))) + y_int_size/5.

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.boxplot([x_fractions, y_fractions, z_fractions], notch=True, patch_artist=True, labels=labels, medianprops=medianprops)
plt.ylabel("Normalized fractional relocalization", fontsize=11)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=8)
plt.savefig("{}_change_by_axis".format(prefix))

sys.exit(0)

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.boxplot([x_length_fractions, y_length_fractions, z_length_fractions], notch=True, patch_artist=True, labels=labels, medianprops=medianprops)
plt.ylabel("Fractional length", fontsize=12)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=0, labelsize=8)
plt.savefig("{}_axis_length".format(prefix))
