import numpy as np
import data_tools as dt
import sys
import os
import linear_algebra as la
import array_tools as at
from scipy import signal as sg
from hmmlearn import hmm
import argparse

def call_peaks(data):
	"""Calls peaks using Gaussian hidden markov model"""
	reshaped_data = data.reshape(-1,1)
	model = hmm.GaussianHMM(n_components=2).fit(reshaped_data)
	scores = model.predict(reshaped_data)

	#determine if peaks are 0 or 1
	zero_indices = np.where(scores == 0)
	one_indices = np.where(scores == 1)
	zero_data = data[zero_indices]
	one_data = data[one_indices]
	if np.mean(zero_data) > np.mean(one_data):
		scores[zero_indices] = 1
		scores[one_indices] = 0

	#find boundaries of peaks
	peaks = []
	in_peak = False
	for i, score in enumerate(scores):
		if in_peak and score == 0:	#end of peak
			in_peak = False
			peak.append(i)
			peaks.append(peak)
		elif not in_peak and score == 1:	#start of peak
			in_peak = True
			peak = [i]

	return peaks

def main():
	parser = argparse.ArgumentParser(description="Identify locus-specific changes between Hi-C datasets")
	parser.add_argument("path1", help="path to intrachromosomal Hi-C BED file 1")
	parser.add_argument("path2", help="path to intrachromosomal Hi-C BED file 2")
	parser.add_argument("-N", default=4, help="number of partitions")
	parser.add_argument("-m", default=0, help="genomic coordinate of centromere")
	parser.add_argument("-s", default=3, help="smoothing parameter for calling relocalization peaks")
	parser.add_argument("-x", default="", help="prefix to minimds.py")
	args = parser.parse_args()

	n = 5

	dir1, name1 = args.path1.split("/")
	dir2, name2 = args.path2.split("/")
	prefix1 = name1.split(".")[0]
	prefix2 = name2.split(".")[0]

	min_error = sys.float_info.max
	for iteration in range(n):
		os.system("python {}minimds.py -m {} -N {} -o {}_ {} {}".format(args.x, args.m, args.N, iteration, args.path1, args.path2))
		
		#load structures
		structure1 = dt.structure_from_file("{}/{}_{}_structure.tsv".format(dir1, iteration, prefix1))	
		structure2 = dt.structure_from_file("{}/{}_{}_structure.tsv".format(dir2, iteration, prefix2))

		#rescale
		structure1.rescale()
		structure2.rescale()

		#make structures compatible
		dt.make_compatible((structure1, structure2))

		#align
		r, t = la.getTransformation(structure1, structure2)
		structure1.transform(r,t)

		#calculate error
		coords1 = np.array(structure1.getCoords())
		coords2 = np.array(structure2.getCoords())
		error = np.mean([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)])
		if error < min_error:
			min_error = error
			best_iteration = iteration

	for iteration in range(n):
		if iteration == best_iteration:
			#load structures
			structure1 = dt.structure_from_file("{}/{}_{}_structure.tsv".format(dir1, iteration, prefix1))	
			structure2 = dt.structure_from_file("{}/{}_{}_structure.tsv".format(dir2, iteration, prefix2))
		else:
			os.system("rm {}/{}_{}_structure.tsv".format(dir1, iteration, prefix1))	
			os.system("rm {}/{}_{}_structure.tsv".format(dir2, iteration, prefix2))		

	#rescale
	structure1.rescale()
	structure2.rescale()

	#make structures compatible
	dt.make_compatible((structure1, structure2))

	#tweak alignment
	r, t = la.getTransformation(structure1, structure2)
	structure1.transform(r,t)

	coords1 = np.array(structure1.getCoords())
	coords2 = np.array(structure2.getCoords())
	dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)]
	print np.mean(dists)

	#smoothed_dists = sg.cwt(dists, sg.ricker, [float(args.s)])[0]
	#dist_peaks = call_peaks(smoothed_dists)
	dist_peaks = sg.find_peaks_cwt(dists, np.arange(1, 20))

	gen_coords = structure1.getGenCoords()

	with open("{}_{}_relocalization.bed".format(prefix1, prefix2), "w") as out:
		for peak in dist_peaks:
			start, end = peak
			peak_dists = dists[start:end]
			max_dist_index = np.argmax(peak_dists) + start
			#out.write("\t".join(("{}".format(structure1.chrom.name), str(gen_coords[start]), str(gen_coords[end]), str(gen_coords[max_dist_index]))))
			out.write("\t".join(("{}".format(structure1.chrom.name), str(gen_coords[max_dist_index]), str(gen_coords[max_dist_index] + structure1.chrom.res)))
			out.write("\n")
		out.close()

if __name__ == "__main__":
	main()
