import data_tools as dt
import numpy as np
import compartment_analysis as ca
from scipy import stats as st
from sklearn import svm
import linear_algebra as la

def calculate_compartment_fraction(structure1, structure2, path1, path2):
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

	x_mean = np.mean(np.abs(x_diffs))/x_length
	y_mean = np.mean(np.abs(y_diffs))/y_length
	z_mean = np.mean(np.abs(z_diffs))/z_length

	return z_mean/(x_mean + y_mean + z_mean)
