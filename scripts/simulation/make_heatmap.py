import sys
sys.path.append("../multimds")
import data_tools as dt
import heatmap as hm

num = sys.argv[1]

path = "sim{}_chr21_100kb.bed".format(num)
struct = dt.structureFromBed(path)
mat = dt.matFromBed(path, struct)

tads = np.loadtxt("sim{}_tads.tsv".format(num))
tad_indices = [(struct.get_rel_index(start), struct.get_rel_index(end)) for start,end in tads]

hm.heatMapFromMat(mat, maxvalue=50, tads=tad_indices)
