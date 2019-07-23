import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
import plotting as plot

cell_type = sys.argv[1]
res_kb = int(sys.argv[2])

struct = dt.structure_from_file("{}_21_{}kb_structure.tsv".format(cell_type, res_kb))

#truncate
start = 45000000
index = struct.chrom.getAbsoluteIndex(start)
struct.points = struct.points[index:len(struct.points)]
struct.chrom.minPos = start
for i in range(len(struct.points)):
	if struct.points[i] != 0:
		struct.points[i].absolute_index -= index
struct.set_rel_indices()

colors = np.zeros_like(struct.getPoints(), dtype=int)
colors[struct.get_rel_index(46900000):struct.get_rel_index(46950000)] = 2
colors[struct.get_rel_index(47475000)] = 1

plot.plot_structure_interactive(struct, colors, colormap="brg")
