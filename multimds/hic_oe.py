import numpy as np

def oe(mat, struct):
	n = len(mat)

	tots = np.zeros(n-1)
	counts = np.zeros(n-1)

	points = struct.getPoints()
	for i in range(n):
		for j in range(i):
			observed = mat[i,j]
			if observed != 0:
				s = points[i].absolute_index - points[j].absolute_index
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
				val = observed/expected
				oe_mat[i,j] = val
				oe_mat[j,i] = val

	return oe_mat

def get_expected(mat):
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

	return avgs
