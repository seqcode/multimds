import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import compartment_analysis as ca
import os
import linear_algebra as la
import array_tools as at

def format_celltype(cell_type):
	formatted = cell_type.split("_")[0]
	return formatted[0].upper() + formatted[1:len(formatted)].lower()

cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
chrom = sys.argv[3]
res = int(sys.argv[4])
res_kb = res/1000

path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

#os.system("python ../multimds.py --full {} {}".format(path1, path2))

#load structures
structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))
structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))
#structure1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))	
#structure2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

dists = np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(structure1.getCoords(), structure2.getCoords())])

#compartments
contacts1 = dt.matFromBed(path1, structure1)
contacts2 = dt.matFromBed(path2, structure2)

active1 = np.array(np.loadtxt("binding_data/{}_{}_{}kb_active_coverage.bed".format(format_celltype(cell_type1), chrom, res_kb), dtype=object)[:,6], dtype=float)
bin_nums1 = structure1.chrom.minPos/structure1.chrom.res + structure1.nonzero_abs_indices()
active1 = active1[bin_nums1]
active2 = np.array(np.loadtxt("binding_data/{}_{}_{}kb_active_coverage.bed".format(format_celltype(cell_type2), chrom, res_kb), dtype=object)[:,6], dtype=float)
bin_nums2 = structure2.chrom.minPos/structure2.chrom.res + structure2.nonzero_abs_indices()
active2 = active2[bin_nums2]

compartments1 = np.array(ca.get_compartments(contacts1, active1))
compartments2 = np.array(ca.get_compartments(contacts2, active2))

n = 50	#number of top peaks

indices = np.flip(np.argsort(dists), 0)
sorted_dists = dists[indices][0:n]
sorted_comps1 = compartments1[indices][0:n]
sorted_comps2 = compartments2[indices][0:n]
gen_coords = np.array(structure1.getGenCoords())[indices][0:n]

with open("{}_A_relocalization.bed".format(chrom), "w") as out:
	for gen_coord, dist, comp1, comp2 in zip(gen_coords, sorted_dists, sorted_comps1, sorted_comps2):
		if np.abs(comp1 - comp2) < 0.2 and comp1 > 0 and comp2 > 0:
			out.write("\t".join((structure1.chrom.name, str(gen_coord), str(gen_coord + structure1.chrom.res))))
			out.write("\n")
	out.close()
