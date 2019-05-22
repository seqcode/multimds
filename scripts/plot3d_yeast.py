import sys
sys.path.append("..")
import data_tools as dt
import os
import numpy as np
from mayavi import mlab

gene_name = sys.argv[1]
chrom_num = sys.argv[2]
gene_loc = int(sys.argv[3])
strain = sys.argv[4]
res_kb = 32

chrom_name = "{}_{}".format(strain, chrom_num)
os.system("python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_{}_{}kb.bed hic_data/galactose_{}_{}kb.bed".format(chrom_name, res_kb, chrom_name, res_kb))
struct1 = dt.structure_from_file("ctrl_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
struct2 = dt.structure_from_file("galactose_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
coords1 = np.array(struct1.getCoords())
coords2 = np.array(struct2.getCoords())

colors = np.zeros_like(struct1.getPoints(), dtype=int)
colors[struct1.get_rel_index(gene_loc)] = 1

mlab.figure(bgcolor=(1,1,1))
line = mlab.plot3d(coords1[:,0], coords1[:,1], coords1[:,2], colors)
lut = line.module_manager.scalar_lut_manager.lut.table.to_array()
lut[0] = (0,0,255,128)
lut[1:len(lut)] = (255,0,0,128)
line.module_manager.scalar_lut_manager.lut.table = lut
line = mlab.plot3d(coords2[:,0], coords2[:,1], coords2[:,2], colors)
lut[0] = (0,255,0,128)
lut[1:len(lut)] = (255,0,0,128)
line.module_manager.scalar_lut_manager.lut.table = lut
mlab.show()
