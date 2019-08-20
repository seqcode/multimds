import sys
from multimds import multimds

celltype1 = sys.argv[1]
celltype2 = sys.argv[2]

multimds.full_mds("hic_data/{}_21_100kb.bed".format(celltype1), "hic_data/{}_21_100kb.bed".format(celltype2))
