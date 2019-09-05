from sklearn.decomposition import PCA
import numpy as np
from scipy import stats as st
from .data_tools import *
from sklearn import svm
from .linear_algebra import *
from .hic_oe import oe

def cor(mat):
	"""Correlation of rows with columns of mat"""
	n = len(mat)
	cor_mat = np.zeros_like(mat)

	for i in range(n):
		for j in range(i):
			r, p = st.pearsonr(mat[i], mat[j])
			cor_mat[i,j] = r
			cor_mat[j,i] = r

	return cor_mat

def get_compartments(mat, enrichments=None, active=True):
	"""From Lieberman-Aiden et al (2009)"""
	oe_mat = oe(mat)
	cor_mat = cor(oe_mat)
	pca = PCA(n_components=1)
	pca.fit(cor_mat)
	scores = pca.fit_transform(cor_mat)[:,0]

	#enforce positive score = active chromatin
	if enrichments is not None:
		r, p = st.pearsonr(scores, enrichments)
		if active and r < 0:
			scores = -scores
		elif not active and r > 0:
			scores = -scores
	
	#normalize
	max_val = max(scores)
	min_val = -min(scores)
	for i, score in enumerate(scores):
		if score > 0:
			scores[i] = score/max_val
		else:
			scores[i] = score/min_val
	
	return scores

def load_enrichments(path, structure, column):
	enrichments = np.array(np.loadtxt(path, dtype=object)[:,column], dtype=float)
	bin_nums = structure.nonzero_abs_indices() + structure.chrom.minPos/structure.chrom.res
	return enrichments[bin_nums]

def calculate_compartment_fraction(structure1, structure2, path1, path2, size1=None, size2=None):
	#compartments
	contacts1 = matFromBed(path1, size1, structure1)
	contacts2 = matFromBed(path2, size2, structure2)

	compartments1 = np.array(get_compartments(contacts1))
	compartments2 = np.array(get_compartments(contacts2))

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

	transformed_coords1 = np.array(change_coordinate_system(coef, coords1))
	transformed_coords2 = np.array(change_coordinate_system(coef, coords2))

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
