import sys
sys.path.append("..")
import compartment_analysis as ca
import data_tools as dt
import numpy as np

res = int(sys.argv[1])
res_kb = int(res/1000)

for chrom in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21):
	path1 = "hic_data/GM12878_combined_{}_100kb.bed".format(chrom)
	structure1 = dt.structureFromBed(path1)

	path2 = "hic_data/K562_{}_100kb.bed".format(chrom)
	structure2 = dt.structureFromBed(path2)

	dt.make_compatible((structure1, structure2))

	contacts = dt.matFromBed(path1, structure1)
	enrichments = np.array(np.loadtxt("binding_data/GM12878_{}_100kb_active_coverage.bed".format(chrom), dtype=object)[:,6], dtype=float)
	bin_nums = structure1.nonzero_abs_indices() + structure1.chrom.minPos/structure1.chrom.res
	enrichments = enrichments[bin_nums]
	compartments1 = np.array(ca.get_compartments(contacts, enrichments))
	
	contacts = dt.matFromBed(path2, structure2)
	enrichments = np.array(np.loadtxt("binding_data/K562_{}_100kb_active_coverage.bed".format(chrom), dtype=object)[:,6], dtype=float)
	bin_nums = structure2.nonzero_abs_indices() + structure2.chrom.minPos/structure2.chrom.res
	enrichments = enrichments[bin_nums]
	compartments2 = np.array(ca.get_compartments(contacts, enrichments))

	gen_coords = structure1.getGenCoords()

	with open("A_compartment_{}kb.bed".format(res_kb), "a") as out:
		for gen_coord, compartment1, compartment2 in zip(gen_coords, compartments1, compartments2):
			if compartment1 > 0 and compartment2 > 0 and np.abs(compartment1 - compartment2) < 0.1:
				for i in range(int(100/res_kb)):
					out.write("\t".join((structure1.chrom.name, str(gen_coord + i*res), str(gen_coord + (i+1)*res), str(compartment1), str(compartment2))))
					out.write("\n")
		out.close()
