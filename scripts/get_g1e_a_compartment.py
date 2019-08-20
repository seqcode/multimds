import sys
sys.path.append("..")
import compartment_analysis as ca
import data_tools as dt
import os
import numpy as np

res = int(sys.argv[1])
res_kb = int(res/1000)

os.system("rm g1e_A_compartment_{}kb.bed".format(res_kb))

for chrom in range(1,20):
	path1 = "hic_data/WT-G1E_{}_100kb.bed".format(chrom)
	structure1 = dt.structureFromBed(path1)

	path2 = "hic_data/KO-rep1-G1E_{}_100kb.bed".format(chrom)
	structure2 = dt.structureFromBed(path2)

	dt.make_compatible((structure1, structure2))

	contacts = dt.matFromBed(path1, structure1)
	enrichments = np.array(np.loadtxt("binding_data/WT-G1E_100kb_CTCF_coverage_{}.bed".format(chrom), dtype=object)[:,6], dtype=float)
	bin_nums = structure1.nonzero_abs_indices() + structure1.chrom.minPos/structure1.chrom.res
	enrichments = enrichments[bin_nums]
	compartments1 = np.array(ca.get_compartments(contacts, enrichments))
	
	contacts = dt.matFromBed(path2, structure2)
	enrichments = np.array(np.loadtxt("binding_data/KO-rep1-G1E_100kb_CTCF_coverage_{}.bed".format(chrom), dtype=object)[:,6], dtype=float)
	bin_nums = structure2.nonzero_abs_indices() + structure2.chrom.minPos/structure2.chrom.res
	enrichments = enrichments[bin_nums]
	compartments2 = np.array(ca.get_compartments(contacts, enrichments))

	gen_coords = structure1.getGenCoords()

	with open("g1e_A_compartment_{}kb.bed".format(res_kb), "a") as out:
		for gen_coord, compartment1, compartment2 in zip(gen_coords, compartments1, compartments2):
			if compartment1 > 0 and compartment2 > 0 and np.abs(compartment1 - compartment2) < 0.1:
				for i in range(int(100/res_kb)):
					out.write("\t".join((structure1.chrom.name, str(gen_coord + i*res), str(gen_coord + (i+1)*res))))
					out.write("\n")
		out.close()