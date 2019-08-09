import sys
import numpy as np
from joint_mds import Joint_MDS
import pymp
import multiprocessing as mp
import data_tools as dt
import array_tools as at
import linear_algebra as la
import tools
import tad
from hic_oe import get_expected
import os
from compartment_fraction import calculate_compartment_fraction

def distmat(contactMat, structure, alpha, weight, num_threads):
	assert len(structure.nonzero_abs_indices()) == len(contactMat)

	expected = get_expected(contactMat)
	for i in range(len(contactMat)):
		for j in range(i):
			corrected = (1-weight)*contactMat[i,j] + weight*expected[i-j-1]
			contactMat[i,j] = corrected
			contactMat[j,i] = corrected

	rowsums = np.array([sum(row) for row in contactMat])
	assert len(np.where(rowsums == 0)[0]) == 0 

	distMat = at.contactToDist(contactMat, alpha)

	distMat = distMat/np.mean(distMat)	#normalize
	
	return distMat

def infer_structures(contactMat1, structure1, contactMat2, structure2, alpha, penalty, num_threads, weight):
	"""Infers 3D coordinates for one structure"""
	distMat1 = distmat(contactMat1, structure1, alpha, weight, num_threads)	
	distMat2 = distmat(contactMat2, structure2, alpha, weight, num_threads)

	coords1, coords2 = Joint_MDS(p=penalty, n_components=3, metric=True, random_state1=np.random.RandomState(), random_state2=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=num_threads).fit_transform(distMat1, distMat2)

	structure1.setCoords(coords1)
	structure2.setCoords(coords2)

def full_mds(path1, path2, alpha=4, penalty=0.05, num_threads=3, weight=0.05, prefix=""):
	"""MDS without partitioning"""
	structure1 = dt.structureFromBed(path1)
	structure2 = dt.structureFromBed(path2)
	dt.make_compatible((structure1, structure2))
	contactMat1 = dt.matFromBed(path1, structure1)
	contactMat2 = dt.matFromBed(path2, structure2)
	infer_structures(contactMat1, structure1, contactMat2, structure2, alpha, penalty, num_threads, weight)

	prefix1 = os.path.splitext(os.path.basename(path1))[0]
	structure1.write("{}{}_structure.tsv".format(prefix, prefix1))
	prefix2 = os.path.splitext(os.path.basename(path2))[0]
	structure2.write("{}{}_structure.tsv".format(prefix, prefix2))

	dists = la.calculate_distances(structure1, structure2)
	with open("{}{}_{}_relocalization.bed".format(prefix, prefix1, prefix2), "w") as out:
		for gen_coord, dist in zip(structure1.getGenCoords(), dists):
			out.write("\t".join((structure1.chrom.name, str(gen_coord), str(gen_coord + structure1.chrom.res), str(dist))))
			out.write("\n")
		out.close()

	print("Fractional compartment change: ")
	print(calculate_compartment_fraction(structure1, structure2, path1, path2))

	return structure1, structure2

def create_low_res_structure(path, res_ratio):
	low_chrom = dt.chromFromBed(path)
	low_chrom.res *= res_ratio
	low_chrom.minPos = int(np.floor(float(low_chrom.minPos)/low_chrom.res)) * low_chrom.res	#round
	low_chrom.maxPos = int(np.ceil(float(low_chrom.maxPos)/low_chrom.res)) * low_chrom.res
	return dt.structureFromBed(path, low_chrom)

def create_high_res_structure(path, lowstructure):
	size, res = dt.basicParamsFromBed(path)
	highChrom = dt.ChromParameters(lowstructure.chrom.minPos, lowstructure.chrom.maxPos, res, lowstructure.chrom.name, size)
	return dt.Structure([], [], highChrom, 0)
	
def transform(trueLow, highSubstructure, res_ratio):
	#approximate as low resolution
	inferredLow = dt.highToLow(highSubstructure, res_ratio)

	scaling_factor = la.radius_of_gyration(trueLow)/la.radius_of_gyration(inferredLow)
	for i, point in enumerate(inferredLow.points):
		if point != 0:
			x, y, z = point.pos
			inferredLow.points[i].pos = (x*scaling_factor, y*scaling_factor, z*scaling_factor)
	
	#recover the transformation for inferred from true low structure
	r, t = la.getTransformation(inferredLow, trueLow)
	t /= scaling_factor

	#transform high structure
	highSubstructure.transform(r, t)

def partitioned_mds(path1, path2, prefix="", centromere=0, num_partitions=4, maxmemory=32000000, num_threads=3, alpha=4, res_ratio=10, penalty=0.05, weight=0.05):
	"""Partitions structure into substructures and performs MDS"""
	#create low-res structures
	lowstructure1 = create_low_res_structure(path1, res_ratio)
	lowstructure2 = create_low_res_structure(path2, res_ratio)
	dt.make_compatible((lowstructure1, lowstructure2))

	#get partitions
	n = len(lowstructure1.getPoints())
	if centromere == 0:
		midpoint = int(n/2)
	else:	
		midpoint = lowstructure1.chrom.getAbsoluteIndex(centromere)
	
	assert num_partitions%2 == 0

	partition_size1 = int(np.ceil(float(midpoint)/(num_partitions/2)))
	partition_size2 = int(np.ceil(float(n-midpoint)/(num_partitions/2)))

	lowpartitions = []	#low substructures, defined on absolute indices not relative indices

	for i in range(int(num_partitions/2)):
		lowpartitions.append((i*partition_size1, min(((i+1)*partition_size1), midpoint)))

	for i in range(int(num_partitions/2)):
		lowpartitions.append((midpoint + i*partition_size2, min((midpoint + (i+1)*partition_size2), n-1)))

	lowpartitions = np.array(lowpartitions)

	low_contactMat1 = dt.matFromBed(path1, lowstructure1)
	low_contactMat2 = dt.matFromBed(path2, lowstructure2)

	tad.substructuresFromAbsoluteTads(lowstructure1, lowpartitions)
	tad.substructuresFromAbsoluteTads(lowstructure2, lowpartitions)

	#create high-res chroms
	size1, res1 = dt.basicParamsFromBed(path1)
	highChrom1 = dt.ChromParameters(lowstructure1.chrom.minPos, lowstructure1.chrom.maxPos, res1, lowstructure1.chrom.name, size1)
	size2, res2 = dt.basicParamsFromBed(path2)
	highChrom2 = dt.ChromParameters(lowstructure2.chrom.minPos, lowstructure2.chrom.maxPos, res2, lowstructure2.chrom.name, size2)

	#initialize high-res substructures
	high_substructures1 = []
	high_substructures2 = []
	low_gen_coords = lowstructure1.getGenCoords()
	offset1 = 0 #initialize
	offset2 = 0
	for partition in lowpartitions:
		start_gen_coord = low_gen_coords[partition[0]]
		end_gen_coord = low_gen_coords[partition[1]]
		high_substructure1 = dt.structureFromBed(path1, highChrom1, start_gen_coord, end_gen_coord, offset1)
		high_substructure2 = dt.structureFromBed(path2, highChrom2, start_gen_coord, end_gen_coord, offset2)
		high_substructures1.append(high_substructure1)
		high_substructures2.append(high_substructure2)
		offset1 += (len(high_substructure1.points) - 1)	#update
		offset2 += (len(high_substructure2.points) - 1)	#update
	
	for high_substructure1, high_substructure2 in zip(high_substructures1, high_substructures2):
		dt.make_points_compatible((high_substructure1, high_substructure2))

	highstructure1 = dt.Structure([], high_substructures1, highChrom1, 0)
	highstructure2 = dt.Structure([], high_substructures2, highChrom2, 0)

	infer_structures(low_contactMat1, lowstructure1, low_contactMat2, lowstructure2, alpha, penalty, num_threads, weight)
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
			trueLow1 = lowSubstructures1[substructurenum]
			trueLow2 = lowSubstructures2[substructurenum]

			#joint MDS
			structure_contactMat1 = dt.matFromBed(path1, highSubstructure1)	#contact matrix for this structure only
			structure_contactMat2 = dt.matFromBed(path2, highSubstructure2)	#contact matrix for this structure only

			infer_structures(structure_contactMat1, highSubstructure1, structure_contactMat2, highSubstructure2, 2.5, penalty, num_threads, weight)

			transform(trueLow1, highSubstructure1, res_ratio)
			transform(trueLow2, highSubstructure2, res_ratio)
	
			highSubstructures1[substructurenum] = highSubstructure1
			highSubstructures2[substructurenum] = highSubstructure2

			print("MDS performed on structure {} of {}".format(substructurenum + 1, numSubstructures))
		
	highstructure1.setstructures(highSubstructures1)
	highstructure2.setstructures(highSubstructures2)

	highstructure1.set_rel_indices()
	highstructure2.set_rel_indices()

	return highstructure1, highstructure2
