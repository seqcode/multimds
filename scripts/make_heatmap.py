import sys
from multimds import data_tools as dt
import heatmap as hm

num = sys.argv[1]

path = "sim{}_chr21_100kb.bed".format(num)
size = dt.size_from_bed(path)
struct = dt.structureFromBed(path, size)
mat = dt.matFromBed(path, size, struct)

tads = np.loadtxt("sim{}_tads.tsv".format(num))
tad_indices = [(struct.get_rel_index(start), struct.get_rel_index(end)) for start,end in tads]

hm.heatMapFromMat(mat, maxvalue=50, tads=tad_indices)
