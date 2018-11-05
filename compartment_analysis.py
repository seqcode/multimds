import array_tools as at
from sklearn.decomposition import PCA
import numpy as np
from scipy import stats as st
import pymp
import sys

def oe(mat, num_threads):
	n = len(mat)

	tots = np.zeros(n-1)
	counts = np.zeros(n-1)
	
	pymp_mat = pymp.shared.array((n,n))
	pymp_mat[0:n,0:n] = mat[0:n,0:n]
	tots = pymp.shared.list(tots)
	counts = pymp.shared.list(counts)

	partition_size = n/num_threads

	with pymp.Parallel(num_threads) as p:
		for thread_num in p.range(num_threads):
			start = partition_size*thread_num
			if thread_num == num_threads-1:
				end = n
			else:
				end = partition_size*(thread_num+1)
			for i in range(start, end):
				for j in range(i):
					observed = pymp_mat[i,j]
					if observed != 0:
						s = i - j
						tots[s - 1] += observed
						counts[s - 1] += 1

	print "binned data"

	avgs = np.zeros(n-1)
	avgs = pymp.shared.list(avgs)
	counts = np.array(counts)
	tots = np.array(tots)

	n = len(counts)
	partition_size = n/num_threads

	with pymp.Parallel(num_threads) as p:
		for thread_num in p.range(num_threads):
                        start = partition_size*thread_num
                        if thread_num == num_threads-1:
                                end = n
                        else:
                                end = partition_size*(thread_num+1)
                        for i in range(start, end):
				if counts[i] != 0:
					avgs[i] = tots[i]/counts[i]

	print "calculated expected"

	n = len(mat)
	partition_size = n/num_threads

	#oe_mat = np.zeros_like(mat)
	oe_mat = pymp.shared.array((n,n))
	avgs = np.array(avgs)

	with pymp.Parallel(num_threads) as p:
                for thread_num in p.range(num_threads):
                        start = partition_size*thread_num
                        if thread_num == num_threads-1:
                                end = n
                        else:
                                end = partition_size*(thread_num+1)
                        for i in range(start, end):
				for j in range(i):
					observed = mat[i,j]
					s = i-j
					expected = avgs[s-1]
					if expected != 0:	
						oe_mat[i,j] = observed/expected

	return oe_mat

def cor(mat, num_threads):
	"""Correlation of rows with columns of mat"""
	n = len(mat)
	#cor_mat = np.zeros_like(mat)
	cor_mat = pymp.shared.array((n,n))
	pymp_mat = pymp.shared.array((n,n))
	pymp_mat[0:n,0:n] = mat[0:n,0:n]

	partition_size = n/num_threads
	with pymp.Parallel(num_threads) as p:
                for thread_num in p.range(num_threads):
                        start = partition_size*thread_num
                        if thread_num == num_threads-1:
                                end = n
                        else:
                                end = partition_size*(thread_num+1)
                        for i in range(start, end):
				for j in range(i):
					r, p = st.pearsonr(pymp_mat[i], pymp_mat[j])
					cor_mat[i,j] = r

	return cor_mat

def get_compartments(mat, num_threads, enrichments=None, active=True):
	"""From Lieberman-Aiden et al (2009)"""
	oe_mat = oe(mat, num_threads)
	print "observed/expected matrix calculated"
	oe_mat = at.makeSymmetric(oe_mat, num_threads)
	print "matrix made symmetric"
	cor_mat = cor(oe_mat, num_threads)
	print "correlation matrix calculated"
	cor_mat = at.makeSymmetric(cor_mat, num_threads)
	print "matrix made symmetric"
	pca = PCA(n_components=1)
	pca.fit(cor_mat)
	print "PCA calculated"
	scores = pca.fit_transform(cor_mat)[:,0]
	print "compartment scores calculated"

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
	print "scores normalized"
	
	return scores

def load_enrichments(path, structure, column):
	enrichments = np.array(np.loadtxt(path, dtype=object)[:,column], dtype=float)
	bin_nums = structure.nonzero_abs_indices() + structure.chrom.minPos/structure.chrom.res
	return enrichments[bin_nums]
