import os
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
#import plotting as plot
from mayavi import mlab

os.system("python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_Suva_13_32kb.bed hic_data/galactose_Suva_13_32kb.bed")
struct1 = dt.structure_from_file("ctrl_Suva_13_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Suva_13_32kb_structure.tsv")

colors = np.zeros_like(struct1.getPoints(), dtype=int)
colors[struct1.get_rel_index(852000)] = 1

#plot.plot_structures_interactive((struct1, struct2), (colors, colors))

mlab.close(all=True)
mlab.figure(bgcolor=(1,1,1))
coords1 = np.array(struct1.getCoords())
mlab.plot3d(coords1[:,0], coords1[:,1], coords1[:,2], colors, colormap="autumn")
coords2 = np.array(struct2.getCoords())
mlab.plot3d(coords2[:,0], coords2[:,1], coords2[:,2], colors, colormap="spring")
mlab.savefig("sup2")
