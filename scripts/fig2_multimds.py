from multimds import multimds
import sys

iteration = sys.argv[1]
chrom = sys.argv[2]
strain = sys.argv[3]

multimds.full_mds("hic_data/ctrl_{}_{}_32kb.bed".format(strain, chrom), "hic_data/galactose_{}_{}_32kb.bed".format(strain, chrom), weight=0, penalty=0.1, output=iteration+"_")
