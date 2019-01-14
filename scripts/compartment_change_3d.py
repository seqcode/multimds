import os
import sys
sys.path.append("..")
import plotting as plot
import data_tools as dt
import numpy as np
import compartment_analysis as ca
from sklearn import svm
import linear_algebra as la

path1 = "hic_data/GM12878_combined_21_100kb.bed"
path2 = "hic_data/K562_21_100kb.bed"
os.system("python ../multimds.py {} {}".format(path1, path2))

struct1 = dt.structure_from_file("GM12878_combined_21_100kb_structure.tsv")	
struct2 = dt.structure_from_file("K562_21_100kb_structure.tsv")

print(struct1.chrom.getGenCoord(325))
print(struct1.chrom.getGenCoord(345))
sys.exit(0)

contacts1 = dt.matFromBed(path1, struct1)
enrichments1 = np.loadtxt("binding_data/Gm12878_21_100kb_active_coverage.bed", usecols=6)
bin_nums1 = struct1.nonzero_abs_indices() + int(struct1.chrom.minPos/struct1.chrom.res)
enrichments1 = enrichments1[bin_nums1]
comps1 = np.array(ca.get_compartments(contacts1, enrichments1))

contacts2 = dt.matFromBed(path2, struct2)
enrichments2 = np.loadtxt("binding_data/K562_21_100kb_active_coverage.bed", usecols=6)
bin_nums2 = struct2.nonzero_abs_indices() + int(struct2.chrom.minPos/struct2.chrom.res)
enrichments2 = enrichments2[bin_nums2]
comps2 = np.array(ca.get_compartments(contacts2, enrichments2))

coords1 = struct1.getCoords()
coords2 = struct2.getCoords()
coords = np.concatenate((coords1, coords2))
compartments = np.concatenate((comps1, comps2))
clf = svm.LinearSVR()
clf.fit(coords, compartments)
coef = clf.coef_

transformed_coords1 = np.array(la.change_coordinate_system(coef, coords1))
transformed_coords2 = np.array(la.change_coordinate_system(coef, coords2))
struct1.setCoords(transformed_coords1)
struct2.setCoords(transformed_coords2)

struct1.subsamplePoints(325, 345)
struct2.subsamplePoints(325, 345)

colors = np.zeros_like(struct1.getPoints(), dtype=int)
index1 = struct1.get_rel_index(42700000)
index2 = struct1.get_rel_index(42900000)
colors[index1:index2] = 1

plot.plot_structures_interactive((struct1, struct2), colors=((0,1,0), (0,1,238./255)), radius=0.01, out_path="compartment_change_3d.png")
plot.plot_structures_interactive((struct1, struct2), (colors, colors), radius=0.01, out_path="compartment_change_3d_mx.png")
