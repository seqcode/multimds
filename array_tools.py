import numpy as np
import pymp

def contactToDist(contactMat, alpha):
	"""Convert contact matrix to distance matrix."""
	distMat = np.zeros_like(contactMat)
	numRows = len(contactMat)
	for i in range(numRows):
		for j in range(i+1):
			if contactMat[i,j] != 0:
				distMat[i,j] = contactMat[i,j]**(-1./alpha)
	return distMat

def makeSymmetric(mat, num_threads):
	"""Make below-diagonal matrix symmetric, in place"""
	n = len(mat)
	pymp_mat = pymp.shared.array((n,n))
	pymp_mat[0:n,0:n] = mat[0:n,0:n]	#initialize with data

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
					pymp_mat[j,i] = pymp_mat[i,j]
	return pymp_mat
