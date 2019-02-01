import sys
sys.path.append("..")
import os
import data_tools as dt
import compartment_analysis as ca
import numpy as np
from sklearn import svm
import linear_algebra as la
from mayavi import mlab

path = "hic_data/GM12878_combined_21_100kb.bed"
os.system("python ../minimds.py {}".format(path))
struct = dt.structure_from_file("GM12878_combined_21_100kb_structure.tsv")

new_start = struct.chrom.getAbsoluteIndex(15000000)
struct.subsamplePoints(new_start, len(struct.points)-3)

#compartments
contacts = dt.matFromBed(path, struct)

compartments = np.array(ca.get_compartments(contacts, 1))

#SVR
coords = struct.getCoords()
clf = svm.LinearSVR()
clf.fit(coords, compartments)
coef = clf.coef_

transformed_coords = np.array(la.change_coordinate_system(coef, coords))
xs = transformed_coords[:,0]
min_x = min(xs)
max_x = max(xs)
x_range = max_x - min_x
ys = transformed_coords[:,1]
min_y = min(ys)
max_y = max(ys)
y_range = max_y - min_y
zs = transformed_coords[:,2]
min_z = min(zs)
max_z = max(zs)

mlab.figure(bgcolor=(1,1,1))
mlab.plot3d(xs, ys, zs, compartments, colormap="bwr")
x_coord = max_x + x_range/10
y_coord = max_y + y_range/10
mlab.quiver3d([0], [0], [1], extent=[0, 0, 0, 0, min_z, max_z], color=(0,0,0), line_width=8)
mlab.savefig("compartments_3D.png")
mlab.show()