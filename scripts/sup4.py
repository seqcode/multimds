from sklearn import svm
import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot
import compartment_analysis as ca
import numpy as np
import linear_algebra as la

struct = dt.structure_from_file("hic_data/GM12878_combined_21_100kb_structure.tsv")
mat = dt.matFromBed("hic_data/GM12878_combined_21_100kb.bed", struct)
comps = ca.get_compartments(mat)

coords = struct.getCoords()
clf = svm.LinearSVR()
clf.fit(coords, comps)
coef = clf.coef_
transformed_coords = np.array(la.change_coordinate_system(coef, coords))
struct.setCoords(transformed_coords)
plot.plot_structure_interactive(struct, comps)
