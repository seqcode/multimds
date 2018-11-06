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
#import plotting as plot

res_kb = 100
cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
chroms = range(1, int(sys.argv[3]))

x_means = []
y_means = []
z_means = []
x_lengths = []
y_lengths = []
z_lengths = []

for chrom in chroms:
	path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
	path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)
	
	if os.path.isfile(path1) and os.path.isfile(path2):
		os.system("python ../multimds.py --full -w 0 {} {}".format(path1, path2))
		structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))
		structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

		#plot.plot_structures_interactive((structure1, structure2))

		#compartments
		contacts1 = dt.matFromBed(path1, structure1)
		contacts2 = dt.matFromBed(path2, structure2)

		at.makeSymmetric(contacts1)
		at.makeSymmetric(contacts2)

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

		x_means.append(np.mean(np.abs(x_diffs)))
		y_means.append(np.mean(np.abs(y_diffs)))
		z_means.append(np.mean(np.abs(z_diffs)))

		#axis lengths
		centroid1 = np.mean(transformed_coords1, axis=0)
		centroid2 = np.mean(transformed_coords2, axis=0)
		x_length1 = np.mean([np.abs(coord1[0] - centroid1[0]) for coord1 in transformed_coords1])
		y_length1 = np.mean([np.abs(coord1[1] - centroid1[1]) for coord1 in transformed_coords1])
		z_length1 = np.mean([np.abs(coord1[2] - centroid1[2]) for coord1 in transformed_coords1])
		x_length2 = np.mean([np.abs(coord2[0] - centroid2[0]) for coord2 in transformed_coords2])
		y_length2 = np.mean([np.abs(coord2[1] - centroid2[1]) for coord2 in transformed_coords2])
		z_length2 = np.mean([np.abs(coord2[2] - centroid2[2]) for coord2 in transformed_coords2])

		x_lengths.append(np.mean((x_length1, x_length2)))
		y_lengths.append(np.mean((y_length1, y_length2)))
		z_lengths.append(np.mean((z_length1, z_length2)))

x_fractions = []
y_fractions = []
z_fractions = []
for x_mean, y_mean, z_mean in zip(x_means, y_means, z_means):
	tot = x_mean + y_mean + z_mean
	x_fractions.append(x_mean/tot)
	y_fractions.append(y_mean/tot)
	z_fractions.append(z_mean/tot)

print(np.mean(z_fractions))

x_length_fractions = []
y_length_fractions = []
z_length_fractions = []
for x_length, y_length, z_length in zip(x_lengths, y_lengths, z_lengths):
	tot = x_length + y_length + z_length
	x_length_fractions.append(x_length/tot)
	y_length_fractions.append(y_length/tot)
	z_length_fractions.append(z_length/tot)

print(x_fractions)
print(y_fractions)
print(z_fractions)

ind = np.arange(len(chroms))  # the x locations for the groups
width = 0.2       # the width of the bars

plt.boxplot([x_fractions, y_fractions, z_fractions], labels=["Orthogonal 1", "Orthogonal 2", "Compartment"])
plt.ylabel("Fractional change")
plt.savefig("{}_{}_change_by_axis".format(cell_type1, cell_type2))
#plt.show()
plt.close()

plt.boxplot([x_length_fractions, y_length_fractions, z_length_fractions], labels=["Orthogonal 1", "Orthogonal 2", "Compartment"])
plt.ylabel("Fractional length")
plt.savefig("{}_{}_axis_length".format(cell_type1, cell_type2))
#plt.show()
plt.close()
