import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import compartment_analysis as ca

cell_types = ("IMR90", "HMEC", "HUVEC", "K562")

loc_comps = np.zeros_like(cell_types, dtype=float)
for i, cell_type in enumerate(cell_types):
	path = "hic_data/{}_21_100kb.bed".format(cell_type)
	struct = dt.structureFromBed(path)
	mat = dt.matFromBed(path, struct)
	comps = ca.get_compartments(mat, struct)
	loc_comps[i] = np.abs(comps[struct.get_rel_index(47400000)])

print loc_comps
