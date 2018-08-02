import numpy as np
import sys

chrom = sys.argv[1]
res = 100000

mat = np.loadtxt("{}_relocalization.tsv".format(chrom))

with open("{}_peaks.bed".format(chrom), "w") as out:
	for i, row in enumerate(mat):
		if i == 0:
			prev = 0
		else:
			prev = mat[i-1,1]
		if i == len(mat) - 1:
			next = 0
		else:
			next = mat[i+1,1]
		diff = row[1]
		if diff > prev and diff > next and row[2] > 0 and row[3] > 0:	#local max in A compartment
			out.write("\t".join(("chr{}".format(chrom), str(int(row[0])), str(int(row[0] + res)), str(diff))))
			out.write("\n")
	out.close()
