import sys
import numpy as np
import pymp
import multiprocessing as mp
import os
from .joint_mds import Joint_MDS
from .compartment_analysis import calculate_compartment_fraction
from .data_tools import *
from .linear_algebra import *
from .tad import *

def infer_structures(path1, structure1, path2, structure2, alpha, penalty, num_threads, weight, size1, size2):
	"""Infers 3D coordinates for one structure"""
	distMat1 = distmat(path1, structure1, size1, alpha, weight)	
	distMat2 = distmat(path2, structure2, size2, alpha, weight)

	coords1, coords2 = Joint_MDS(p=penalty, n_components=3, metric=True, random_state1=np.random.RandomState(), random_state2=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=num_threads).fit_transform(distMat1, distMat2)
	
	structure1.setCoords(coords1)
	structure2.setCoords(coords2)

def full_mds(path1, path2, alpha=4, penalty=0.05, num_threads=3, weight=0.05, prefix=""):
	"""MDS without partitioning"""
	#get file sizes
	size1 = size_from_bed(path1)
	size2 = size_from_bed(path2)

	structure1 = structureFromBed(path1, size1)
	structure2 = structureFromBed(path2, size2)
	make_compatible((structure1, structure2))
	infer_structures(path1, structure1, path2, structure2, alpha, penalty, num_threads, weight, size1, size2)

	prefix1 = os.path.splitext(os.path.basename(path1))[0]
	structure1.write("{}{}_structure.tsv".format(prefix, prefix1))
	prefix2 = os.path.splitext(os.path.basename(path2))[0]
	structure2.write("{}{}_structure.tsv".format(prefix, prefix2))

	dists = calculate_distances(structure1, structure2)
	with open("{}{}_{}_relocalization.bed".format(prefix, prefix1, prefix2), "w") as out:
		for gen_coord, dist in zip(structure1.getGenCoords(), dists):
			out.write("\t".join((structure1.chrom.name, str(gen_coord), str(gen_coord + structure1.chrom.res), str(dist))))
			out.write("\n")
		out.close()

	#print("Fractional compartment change: ")
	#print(calculate_compartment_fraction(structure1, structure2, path1, path2, size1, size2))

	return structure1, structure2

def partitioned_mds(path1, path2, prefix="", num_partitions=4, maxmemory=32000000, num_threads=3, alpha=4, res_ratio=10, penalty=0.05, weight=0.05):
	"""Partitions structure into substructures and performs MDS"""
	#get file sizes
	size1 = size_from_bed(path1)
	size2 = size_from_bed(path2)

	#create chromosomes
	highChrom1 = chromFromBed(path1)
	highChrom2 = chromFromBed(path2)
	highChrom = consensus_chrom((highChrom1, highChrom2))
	lowChrom = highChrom.reduceRes(res_ratio)

	#create low-res structures
	lowstructure1 = structureFromBed(path1, size1, lowChrom)
	lowstructure2 = structureFromBed(path2, size2, lowChrom)
	make_compatible((lowstructure1, lowstructure2))

	#get partitions
	n = len(lowstructure1.getPoints())
	assert num_partitions < n/2
	partition_size = int(n/num_partitions)
	lowpartitions = np.array([(i*partition_size, min(((i+1)*partition_size), n-1)) for i in range(num_partitions)])	#low substructures, defined on absolute indices not relative indices

	substructuresFromAbsoluteTads(lowstructure1, lowpartitions)
	substructuresFromAbsoluteTads(lowstructure2, lowpartitions)

	#initialize high-res substructures
	high_substructures1 = []
	high_substructures2 = []
	low_gen_coords = lowstructure1.getGenCoords()
	offset1 = 0 #initialize
	offset2 = 0
	for partition in lowpartitions:
		start_gen_coord = low_gen_coords[partition[0]]
		end_gen_coord = low_gen_coords[partition[1]]
		high_substructure1 = structureFromBed(path1, size1, highChrom, start_gen_coord, end_gen_coord, offset1)
		high_substructure2 = structureFromBed(path2, size2, highChrom, start_gen_coord, end_gen_coord, offset2)
		high_substructures1.append(high_substructure1)
		high_substructures2.append(high_substructure2)
		offset1 += (len(high_substructure1.points) - 1)	#update
		offset2 += (len(high_substructure2.points) - 1)	#update
	
	for high_substructure1, high_substructure2 in zip(high_substructures1, high_substructures2):
		make_points_compatible((high_substructure1, high_substructure2))

	highstructure1 = Structure([], high_substructures1, highChrom1, 0)
	highstructure2 = Structure([], high_substructures2, highChrom2, 0)

	infer_structures(path1, lowstructure1, path2, lowstructure2, alpha, penalty, num_threads, weight, size1, size2)
	print("Low-resolution MDS complete")

	highSubstructures1 = pymp.shared.list(highstructure1.structures)
	highSubstructures2 = pymp.shared.list(highstructure2.structures)
	lowSubstructures1 = pymp.shared.list(lowstructure1.structures)
	lowSubstructures2 = pymp.shared.list(lowstructure2.structures)

	numSubstructures = len(highstructure1.structures)
	num_threads = min((num_threads, mp.cpu_count(), numSubstructures))	#don't exceed number of requested threads, available threads, or structures
	with pymp.Parallel(num_threads) as p:
		for substructurenum in p.range(numSubstructures):
			highSubstructure1 = highSubstructures1[substructurenum]
			highSubstructure2 = highSubstructures2[substructurenum]

			#joint MDS
			infer_structures(path1, highSubstructure1, path2, highSubstructure2, 2.5, penalty, num_threads, weight, size1, size2)

			transform(lowSubstructures1[substructurenum], highSubstructure1, res_ratio)
			transform(lowSubstructures2[substructurenum], highSubstructure2, res_ratio)

			highSubstructures1[substructurenum] = highSubstructure1
			highSubstructures2[substructurenum] = highSubstructure2

			print("MDS performed on structure {} of {}".format(substructurenum + 1, numSubstructures))
		
	highstructure1.setstructures(highSubstructures1)
	highstructure2.setstructures(highSubstructures2)

	highstructure1.set_rel_indices()
	highstructure2.set_rel_indices()

	prefix1 = os.path.splitext(os.path.basename(path1))[0]
	highstructure1.write("{}{}_structure.tsv".format(prefix, prefix1))
	prefix2 = os.path.splitext(os.path.basename(path2))[0]
	highstructure2.write("{}{}_structure.tsv".format(prefix, prefix2))

	dists = calculate_distances(highstructure1, highstructure2)
	with open("{}{}_{}_relocalization.bed".format(prefix, prefix1, prefix2), "w") as out:
		for gen_coord, dist in zip(highstructure1.getGenCoords(), dists):
			out.write("\t".join((highstructure1.chrom.name, str(gen_coord), str(gen_coord + highstructure1.chrom.res), str(dist))))
			out.write("\n")
		out.close()

	#print("Fractional compartment change: ")
	#print(calculate_compartment_fraction(highstructure1, highstructure2, path1, path2, size1, size2))	

	return highstructure1, highstructure2
