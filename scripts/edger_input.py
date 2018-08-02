import sys
sys.path.append("..")
import data_tools as dt
import array_tools as at
import numpy as np

def compatible_chroms(paths):
	chroms = [dt.chromFromBed(path) for path in paths]
	all_min_pos = [chrom.minPos for chrom in chroms]
	all_max_pos = [chrom.maxPos for chrom in chroms]
	consensus_min = max(all_min_pos)
	consensus_max = min(all_max_pos)
	for chrom in chroms:
		chrom.minPos = consensus_min
		chrom.maxPos = consensus_max
	return chroms

def fullMatFromBed(path, chrom):	
	"""Converts BED file to matrix"""
	numpoints = (chrom.maxPos - chrom.minPos)/chrom.res + 1
	mat = np.zeros((numpoints, numpoints))	

	with open(path) as infile:
		for line in infile:
			line = line.strip().split()	#line as array of strings
			loc1 = int(line[1])
			loc2 = int(line[4])
			index1 = chrom.getPointNum(loc1)
			index2 = chrom.getPointNum(loc2)
			if index1 > index2:
				row = index1
				col = index2
			else:
				row = index2
				col = index1
			mat[row, col] += float(line[6])
		infile.close()

	at.makeSymmetric(mat)

	return mat

res_kb = 100
cell_types = ("K562", "GM12878_primary", "GM12878_replicate")

for chrom_name in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22):
	paths = ["hic_data/{}_{}_{}kb.bed".format(cell_type, chrom_name, res_kb) for cell_type in cell_types]
	chroms = compatible_chroms(paths)

	mats = [fullMatFromBed(path, chrom) for path, chrom in zip(paths, chroms)]

	sum_mat = np.sum(mats, 0)

	with open("chr{}_edgeR_table.tsv".format(chrom_name), "w") as out:
		out.write("Symbol\t")
		out.write("\t".join(cell_types))	#header
		out.write("\n")
		for i in range(len(sum_mat[0])):
			for j in range(i):
				if sum_mat[i,j] != 0:	#at least one element is non-zero
					loc1 = chrom.minPos + chrom.res * j
					loc2 = chrom.minPos + chrom.res * i
					out.write("chr{}:{},{}\t".format(chrom_name, loc1, loc2))	#identifier
					out.write("\t".join([str(mat[i,j]) for mat in mats]))
					out.write("\n")
		out.close() 		
