import sys
sys.path.append("..")
import data_tools as dt
import numpy as np
import linear_algebra as la
from sklearn.manifold import MDS

chrom = sys.argv[1]
res_kb = 100
prefix1 = "GM12878_combined"
prefix2 = "K562"
	
path1 = "hic_data/{}_{}_{}kb.bed".format(prefix1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(prefix2, chrom, res_kb)

structure1 = dt.structureFromBed(path1, None, None)
structure2 = dt.structureFromBed(path2, None, None)

#make structures compatible
dt.make_compatible((structure1, structure2))

#get distance matrices
size1 = dt.size_from_bed(path1)
size2 = dt.size_from_bed(path2)
dists1 = dt.distmat(path1, structure1, size1)
dists2 = dt.distmat(path2, structure2, size2)

#MDS
coords1 = MDS(n_components=3, random_state=np.random.RandomState(), dissimilarity="precomputed", n_jobs=-1).fit_transform(dists1)
coords2 = MDS(n_components=3, random_state=np.random.RandomState(), dissimilarity="precomputed", n_jobs=-1).fit_transform(dists2)

#fill structures
structure1.setCoords(coords1)
structure2.setCoords(coords2)

#align
r, t = la.getTransformation(structure1, structure2)
structure1.transform(r, t)
