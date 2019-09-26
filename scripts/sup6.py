import sys
from multimds import data_tools as dt
import heatmap as hm
import numpy as np

num = sys.argv[1]

path = "sim{}_chr21_100kb.bed".format(num)
size = dt.size_from_bed(path)
struct = dt.structureFromBed(path, size)
mat = dt.matFromBed(path, size, struct)

tad_boundaries = np.loadtxt("subseqs{}.tsv".format(num))
tads = []
for i in range(1, len(tad_boundaries)):
	bound1 = tad_boundaries[i-1]
	bound2 = tad_boundaries[i]
	if bound1 == 0:
		start = 0
	else:
		start = struct.get_rel_index(bound1)
	end = struct.get_rel_index(bound2)	
	tads.append((start,end))

#hm.heatMapFromMat(mat, maxvalue=50, tads=tad_indices, outpath="sup6_{}".format(num))
hm.heatMapFromMat(mat, maxvalue=5, tads=tads, outpath="sup6_{}.svg".format(num))
