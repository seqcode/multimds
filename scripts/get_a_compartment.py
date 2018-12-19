import sys
sys.path.append("..")
import compartment_analysis as ca
import data_tools as dt
import array_tools as at
import os
import numpy as np

res = int(sys.argv[1])
res_kb = res/1000

if os.path.isfile("A_compartment_{}kb.bed".format(res_kb)):
	os.system("rm A_compartment_{}kb.bed".format(res_kb))

for chrom in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22):
	path = "hic_data/GM12878_combined_{}_100kb.bed".format(chrom)
	structure = dt.structureFromBed(path)
	contacts = dt.matFromBed(path, structure)
	at.makeSymmetric(contacts)
	enrichments = np.array(np.loadtxt("binding_data/Gm12878_{}_100kb_active_coverage.bed".format(chrom), dtype=object)[:,6], dtype=float)
	bin_nums = structure.nonzero_abs_indices() + structure.chrom.minPos/structure.chrom.res
	enrichments = enrichments[bin_nums]
	compartments = np.array(ca.get_compartments(contacts, enrichments))
	gen_coords = np.array(structure.getGenCoords())
	a_gen_coords = gen_coords[np.where(compartments > 0)]
	with open("A_compartment_{}kb.bed".format(res_kb), "a") as out:
		for a_gen_coord in a_gen_coords:
			for i in range(structure.chrom.res/res):
				out.write("\t".join((structure.chrom.name, str(a_gen_coord + i*res), str(a_gen_coord + (i+1)*res))))
				out.write("\n")
		out.close()
