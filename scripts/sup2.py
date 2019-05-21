import os
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
from mayavi import mlab

strain = sys.argv[1]
chrom = sys.argv[2]
gene_loc = int(sys.argv[3])

os.system("python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_{}_{}_32kb.bed hic_data/galactose_{}_{}_32kb.bed".format(strain, chrom, strain, chrom))
struct1 = dt.structure_from_file("ctrl_{}_{}_32kb_structure.tsv".format(strain, chrom))
struct2 = dt.structure_from_file("galactose_{}_{}_32kb_structure.tsv".format(strain, chrom))

colors = np.zeros_like(struct1.getPoints(), dtype=int)
colors[struct1.get_rel_index(gene_loc)] = -1

mlab.close(all=True)
mlab.figure(bgcolor=(1,1,1))
coords1 = np.array(struct1.getCoords())
mlab.plot3d(coords1[:,0], coords1[:,1], coords1[:,2], colors, colormap="RdYlBu")
coords2 = np.array(struct2.getCoords())
mlab.plot3d(coords2[:,0], coords2[:,1], coords2[:,2], colors, colormap="RdYlGn")
mlab.show()
