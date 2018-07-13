import sys
import data_tools as dt

res_kb = 100
chrom = sys.argv[1]
cell_type1 = "GM12878_combined"
cell_type2 = "K562"
	
path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

structure1 = dt.structureFromBed(path1)
structure2 = dt.structureFromBed(path2)

dt.make_compatible((structure1, structure2))

print "size\t" + str(len(structure1.getPoints()))
