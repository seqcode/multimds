import array_tools as at
from sklearn.decomposition import PCA
import numpy as np
from scipy import stats as st

def oe(mat):
	n = len(mat)

	tots = np.zeros(n-1)
	counts = np.zeros(n-1)

	for i in range(n):
		for j in range(i):
			observed = mat[i,j]
			if observed != 0:
				s = i - j
				tots[s - 1] += observed
				counts[s - 1] += 1

	avgs = np.zeros(n-1)
	for i, count in enumerate(counts):
		if count != 0:
			avgs[i] = tots[i]/count

	oe_mat = np.zeros_like(mat)

	for i in range(n):
		for j in range(i):
			observed = mat[i,j]
			s = i - j
			expected = avgs[s - 1]
			if expected != 0:	
				oe_mat[i,j] = observed/expected

	return oe_mat

def cor(mat):
	"""Correlation of rows with columns of mat"""
	cor_mat = np.zeros_like(mat)
	for i in range(len(mat)):
		for j in range(i):
			r, p = st.pearsonr(mat[i], mat[j])
			cor_mat[i,j] = r
	return cor_mat

def get_compartments(mat, structure, path, active):
	"""From Lieberman-Aiden et al (2009)"""
	oe_mat = oe(mat)
	at.makeSymmetric(oe_mat)
	cor_mat = cor(oe_mat)
	at.makeSymmetric(cor_mat)
	pca = PCA(n_components=1)
	pca.fit(cor_mat)
	scores = pca.fit_transform(cor_mat)[:,0]

	#enforce positive score = active chromatin
	enrichments = np.array(np.loadtxt(path, dtype=object)[:,6], dtype=float)
	bin_nums = structure.getPointNums() + structure.chrom.minPos/structure.chrom.res
	enrichments = enrichments[bin_nums]
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

def infer_compartments(mat, structure, cell_type, chrom, res_kb):
	"""Default method for my file paths"""
	#re-format cell type
	formatted = cell_type.split("_")[0]
	formatted = formatted[0].upper() + formatted[1:len(formatted)].lower()

	return get_compartments(mat, structure, "binding_data/{}_{}_{}kb_active_coverage.bed".format(formatted, chrom, res_kb), True)
