import sys
sys.path.append("..")
import data_tools as dt
import linear_algebra as la
import numpy as np
from mayavi import mlab

struct1 = dt.structure_from_file("hic_data/GM12878_combined_21_10kb_structure.tsv")
struct2 = dt.structure_from_file("hic_data/K562_21_10kb_structure.tsv")

start = 850
end = 950

coords1 = np.array(struct1.getCoords())
coords2 = np.array(struct2.getCoords())
coords1 = coords1[start:end]
coords2 = coords2[start:end]
gen_coords = struct1.getGenCoords()
dists = np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)])
index1 = np.argmax(dists)
index2 = np.argmin(dists)

mlab.close(all=True)
mlab.figure(bgcolor=(1,1,1))
mlab.plot3d(coords1[:,0], coords1[:,1], coords1[:,2], color=(0.5,0,0), tube_radius=0.005)
mlab.plot3d(coords2[:,0], coords2[:,1], coords2[:,2], color=(0,0.5,0), tube_radius=0.005)
mlab.plot3d((coords1[index1,0], coords2[index1,0]), (coords1[index1,1], coords2[index1,1]), (coords1[index1,2], coords2[index1,2]), color=(1,0,1), tube_radius=0.005)
mlab.plot3d((coords1[index2,0], coords2[index2,0]), (coords1[index2,1], coords2[index2,1]), (coords1[index2,2], coords2[index2,2]), color=(0,0,1), tube_radius=0.005)
mlab.show()
