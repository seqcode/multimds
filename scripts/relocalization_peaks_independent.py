import numpy as np
import sys
from multimds import data_tools as dt
from multimds import compartment_analysis as ca
from multimds import linear_algebra as la
from scipy import signal as sg

def format_celltype(cell_type):
	#formatted = cell_type.split("_")[0]
	#return formatted[0].upper() + formatted[1:len(formatted)].lower()
	return cell_type.split("_")[0]

cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
chrom = sys.argv[3]
res = int(sys.argv[4])
res_kb = res/1000

path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

#load structures
structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_independent_structure.tsv".format(cell_type1, chrom, res_kb))	
structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_independent_structure.tsv".format(cell_type2, chrom, res_kb))

#rescale
structure1.rescale()
structure2.rescale()

#make structures compatible
dt.make_compatible((structure1, structure2))

#align
r,t = la.getTransformation(structure1, structure2)
structure1.transform(r,t)

#calculate changes
coords1 = np.array(structure1.getCoords())
coords2 = np.array(structure2.getCoords())
dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)]

#compartments
chrom1 = dt.chromFromBed(path1)
chrom2 = dt.chromFromBed(path2)
chrom1.res = 100000	#reduce res to reduce RAM usage in compartment calculation
chrom2.res = 100000
chrom1.minPos = int(np.floor(float(chrom1.minPos)/chrom1.res)) * chrom1.res	#round
chrom1.maxPos = int(np.ceil(float(chrom1.maxPos)/chrom1.res)) * chrom1.res
chrom2.minPos = int(np.floor(float(chrom2.minPos)/chrom2.res)) * chrom2.res	#round
chrom2.maxPos = int(np.ceil(float(chrom2.maxPos)/chrom2.res)) * chrom2.res

low_struct1 = dt.structureFromBed(path1, chrom1)
low_struct2 = dt.structureFromBed(path2, chrom2)
dt.make_compatible((low_struct1, low_struct2))
contacts1 = dt.matFromBed(path1, low_struct1)		
contacts2 = dt.matFromBed(path2, low_struct2)

enrichments = np.loadtxt("binding_data/{}_{}_100kb_active_coverage.bed".format(format_celltype(cell_type1), chrom), usecols=6)
bin_nums = low_struct1.nonzero_abs_indices() + low_struct1.chrom.minPos/low_struct1.chrom.res
enrichments = enrichments[bin_nums]
compartments1 = np.array(ca.get_compartments(contacts1, enrichments))

enrichments = np.loadtxt("binding_data/{}_{}_100kb_active_coverage.bed".format(format_celltype(cell_type2), chrom), usecols=6)
bin_nums = low_struct2.nonzero_abs_indices() + low_struct2.chrom.minPos/low_struct2.chrom.res
enrichments = enrichments[bin_nums]
compartments2 = np.array(ca.get_compartments(contacts2, enrichments))

gen_coords = structure1.getGenCoords()

compartment_diffs = np.abs(compartments1 - compartments2)

dist_peaks = sg.find_peaks_cwt(dists, np.arange(1,10))

high_coords = structure1.getGenCoords()
low_coords = low_struct1.getGenCoords()

high_comps1 = []
high_comps2 = []

for high_coord in high_coords:
	low_index1 = low_struct1.get_rel_index(high_coord)
	assert low_index1 is not None
	low_comp1 = compartments1[low_index1]
	high_comps1.append(low_comp1)
	low_index2 = low_struct2.get_rel_index(high_coord)
	assert low_index2 is not None
	low_comp2 = compartments2[low_index2]
	high_comps2.append(low_comp2)

with open("{}_dist_peaks_independent.bed".format(chrom), "w") as out:
	for peak in dist_peaks:
		out.write("\t".join((structure1.chrom.name, str(high_coords[peak]), str(high_coords[peak] + structure1.chrom.res), str(high_comps1[peak]), str(high_comps2[peak]))))
		out.write("\n")
	out.close()
