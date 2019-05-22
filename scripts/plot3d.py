import sys
sys.path.append("..")
import data_tools as dt
import os
import numpy as np
import plotting as plot

cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
res_kb = 25
os.system("python ../multimds.py hic_data/{}_21_{}kb.bed hic_data/{}_21_{}kb.bed".format(cell_type1, res_kb, cell_type2, res_kb))

struct1 = dt.structure_from_file("{}_21_{}kb_structure.tsv".format(cell_type1, res_kb))

#truncate
start = 45000000
index = struct1.chrom.getAbsoluteIndex(start)
struct1.points = struct1.points[index:len(struct1.points)]
struct1.chrom.minPos = start
for i in range(len(struct1.points)):
	if struct1.points[i] != 0:
		struct1.points[i].absolute_index -= index
struct1.set_rel_indices()

colors = np.zeros_like(struct1.getPoints(), dtype=int)
colors[struct1.get_rel_index(46900000):struct1.get_rel_index(46950000)] = 1
colors[struct1.get_rel_index(47475000)] = 2

plot.plot_structure_interactive(struct1, colors)

struct2 = dt.structure_from_file("{}_21_{}kb_structure.tsv".format(cell_type2, res_kb))

#truncate
index = struct2.chrom.getAbsoluteIndex(start)
struct2.points = struct2.points[index:len(struct2.points)]
struct2.chrom.minPos = start
for i in range(len(struct2.points)):
	if struct2.points[i] != 0:
		struct2.points[i].absolute_index -= index
struct2.set_rel_indices()

colors = np.zeros_like(struct2.getPoints(), dtype=int)
colors[struct2.get_rel_index(46900000):struct2.get_rel_index(46950000)] = 1
colors[struct2.get_rel_index(47475000)] = 2

plot.plot_structure_interactive(struct2, colors)
