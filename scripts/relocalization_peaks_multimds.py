from multimds import multimds
import sys

chrom = sys.argv[1]

multimds.full_mds("hic_data/GM12878_combined_{}_100kb.bed".format(chrom), "hic_data/K562_{}_100kb.bed".format(chrom))
