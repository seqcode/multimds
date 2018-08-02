import sys
import numpy as np
from joint_mds import Joint_MDS
import argparse
import pymp
import multiprocessing as mp
import data_tools as dt
import array_tools as at
import linear_algebra as la
import tools
import tad

def distmat(contactMat, structure, alpha):
	assert len(structure.nonzero_abs_indices()) == len(contactMat)

	at.makeSymmetric(contactMat)
	rowsums = np.array([sum(row) for row in contactMat])
	assert len(np.where(rowsums == 0)[0]) == 0 

	distMat = at.contactToDist(contactMat, alpha)
	at.makeSymmetric(distMat)

	distMat = distMat/np.mean(distMat)	#normalize
	
	return distMat

def infer_structures(contactMat1, structure1, contactMat2, structure2, alpha, penalty, num_threads):
	"""Infers 3D coordinates for one structure"""
	distMat1 = distmat(contactMat1, structure1, alpha)	
	distMat2 = distmat(contactMat2, structure2, alpha)

	coords1, coords2 = Joint_MDS(p=penalty, n_components=3, metric=True, random_state1=np.random.RandomState(), random_state2=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=num_threads).fit_transform(distMat1, distMat2)

	structure1.setCoords(coords1)
	structure2.setCoords(coords2)

def fullMDS(path1, path2, alpha, penalty, num_threads):
	"""MDS without partitioning"""
	structure1 = dt.structureFromBed(path1)
	structure2 = dt.structureFromBed(path2)
	dt.make_compatible((structure1, structure2))
	contactMat1 = dt.matFromBed(path1, structure1)
	contactMat2 = dt.matFromBed(path2, structure2)
	infer_structures(contactMat1, structure1, contactMat2, structure2, alpha, penalty, num_threads)
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

def partitionedMDS(path1, path2, args):
	"""Partitions structure into substructures and performs MDS"""
	centromere = args[0]
	num_partitions = args[1]
	maxmemory = args[2]
	num_threads = args[3]
	alpha = args[4]
	res_ratio = args[5]
	penalty = args[6]

	#create low-res structures
	lowstructure1 = create_low_res_structure(path1, res_ratio)
	lowstructure2 = create_low_res_structure(path2, res_ratio)
	dt.make_compatible((lowstructure1, lowstructure2))

	#get partitions
	n = len(lowstructure1.getPoints())
	if centromere == 0:
		midpoint = n/2
	else:	
		midpoint = lowstructure1.chrom.getAbsoluteIndex(centromere)
	
	assert num_partitions%2 == 0

	partition_size1 = int(np.ceil(float(midpoint)/(num_partitions/2)))
	partition_size2 = int(np.ceil(float(n-midpoint)/(num_partitions/2)))

	lowpartitions = []	#low substructures, defined on absolute indices not relative indices

	for i in range(num_partitions/2):
		lowpartitions.append((i*partition_size1, min(((i+1)*partition_size1), midpoint)))

	for i in range(num_partitions/2):
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

	infer_structures(low_contactMat1, lowstructure1, low_contactMat2, lowstructure2, alpha, penalty, num_threads)
	print "Low-resolution MDS complete"

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

			infer_structures(structure_contactMat1, highSubstructure1, structure_contactMat2, highSubstructure2, 2.5, penalty, num_threads)

			transform(trueLow1, highSubstructure1, res_ratio)
			transform(trueLow2, highSubstructure2, res_ratio)
	
			highSubstructures1[substructurenum] = highSubstructure1
			highSubstructures2[substructurenum] = highSubstructure2

			print "MDS performed on structure {} of {}".format(substructurenum + 1, numSubstructures)
		
	highstructure1.setstructures(highSubstructures1)
	highstructure2.setstructures(highSubstructures2)

	return highstructure1, highstructure2

def main():
	parser = argparse.ArgumentParser(description="Jointly reconstruct 3D coordinates from two normalized intrachromosomal Hi-C BED files.")
	parser.add_argument("path1", help="path to first intrachromosomal Hi-C BED file")
	parser.add_argument("path2", help="path to second intrachromosomal Hi-C BED file")
	parser.add_argument("--full", action="store_true", help="use full MDS (default: partitioned MDS)")
	parser.add_argument("-l", type=int, help="low resolution/high resolution", default=10)
	parser.add_argument("-o", help="output file prefix")
	parser.add_argument("-r", default=32000000, help="maximum RAM to use (in kb)")
	parser.add_argument("-n", type=int, default=3, help="number of threads")
	parser.add_argument("-a", type=float, default=4, help="alpha factor for converting contact frequencies to physical distances")
	parser.add_argument("-P", type=float, default=0.05, help="joint MDS penalty")
	parser.add_argument("-m", type=int, default=0, help="midpoint (usually centromere) for partitioning")
	parser.add_argument("-N", type=int, default=2, help="number of partitions")
	args = parser.parse_args()

	if args.full:	#not partitioned
		structure1, structure2 = fullMDS(args.path1, args.path2, args.a, args.P, args.n)
	else:	#partitioned
		params = (args.m, args.N, args.r, args.n, args.a, args.l, args.P)
		names = ("Midpoint", "Number of partitions", "Maximum memory", "Number of threads", "Alpha", "Resolution ratio", "Penalty")
		intervals = ((None, None), (1, None), (0, None), (0, None), (1, None), (1, None), (0, None))
		if not tools.args_are_valid(params, names, intervals):
			sys.exit(0)

		structure1, structure2 = partitionedMDS(args.path1, args.path2, params)
	
	if args.o:
		prefix = args.o
	else:
		prefix = ""
	prefix1 = args.path1.split("/")[-1].split(".bed")[0]
	structure1.write("{}/{}{}_structure.tsv".format("/".join(args.path1.split("/")[0:-1]), prefix, prefix1))
	prefix2 = args.path2.split("/")[-1].split(".bed")[0]
	structure2.write("{}/{}{}_structure.tsv".format("/".join(args.path2.split("/")[0:-1]), prefix, prefix2))

if __name__ == "__main__":
	main()
