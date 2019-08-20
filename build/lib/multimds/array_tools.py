import numpy as np

def contactToDist(contactMat, alpha):
	"""Convert contact matrix to distance matrix."""
	distMat = np.zeros_like(contactMat)
	numRows = len(contactMat)
	for i in range(numRows):
		for j in range(i+1):
			if contactMat[i,j] != 0:
				dist = contactMat[i,j]**(-1./alpha)
				distMat[i,j] = dist
				distMat[j,i] = dist
	return distMat
