import os
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot

os.system("python ../multimds.py -P 0.1 -w 0 ctrl_Scer_13_32kb.bed galactose_Scer_13_32kb.bed")
struct1 = dt.structure_from_file("ctrl_Suva_13_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Suva_13_32kb_structure.tsv")

colors = np.zeros_like(struct1.getPoints(), dtype=int)
colors[struct1.get_rel_index(852000)] = 1

plot.plot_structures_interactive((struct1, struct2), (colors, colors))
