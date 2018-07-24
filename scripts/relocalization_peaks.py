import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import compartment_analysis as ca
import os
import linear_algebra as la
import array_tools as at
from scipy import signal as sg
#from hmmlearn import hmm

def normalize(values):
	return np.array(values)/max(values)

def format_celltype(cell_type):
	if cell_type == "KBM7":
		return "K562"	#substitute
	else:
		formatted = cell_type.split("_")[0]
		return formatted[0].upper() + formatted[1:len(formatted)].lower()

#def call_peaks(data):
#	"""Calls peaks using Gaussian hidden markov model"""
#	reshaped_data = data.reshape(-1,1)
#	model = hmm.GaussianHMM(n_components=2).fit(reshaped_data)
#	scores = model.predict(reshaped_data)

	#determine if peaks are 0 or 1
#	zero_indices = np.where(scores == 0)
#	one_indices = np.where(scores == 1)
#	zero_data = data[zero_indices]
#	one_data = data[one_indices]
#	if np.mean(zero_data) > np.mean(one_data):
#		scores[zero_indices] = 1
#		scores[one_indices] = 0

	#find boundaries of peaks
#	peaks = []
#	in_peak = False
#	for i, score in enumerate(scores):
#		if in_peak and score == 0:	#end of peak
#			in_peak = False
#			peak.append(i)
#			peaks.append(peak)
#		elif not in_peak and score == 1:	#start of peak
#			in_peak = True
#			peak = [i]

#	return peaks

res_kb = 100
cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
chrom = sys.argv[3]
centromere = sys.argv[4]
num_partitions = sys.argv[5]
smoothing_parameter = float(sys.argv[6])
n = 1

path1 = "/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}kb_filtered.bed".format(cell_type1, chrom, res_kb)
path2 = "/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}kb_filtered.bed".format(cell_type2, chrom, res_kb)
#path1 = "hic_data/{}_{}_{}kb_filtered.bed".format(cell_type1, chrom, res_kb)
#path2 = "hic_data/{}_{}_{}kb_filtered.bed".format(cell_type2, chrom, res_kb)

min_error = sys.float_info.max
for iteration in range(n):
	os.system("python ../multimds.py -m {} -N {} -o {}_ {} {}".format(centromere, num_partitions, iteration, path1, path2))
		
	#load structures
	structure1 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
	structure2 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type2, chrom, res_kb))

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
		structure1 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
		structure2 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type2, chrom, res_kb))
	else:
		os.system("rm /data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
		os.system("rm /data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type2, chrom, res_kb))	

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
dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)]
print np.mean(dists)

#compartments
contacts1 = dt.matFromBed(path1, structure1)
contacts2 = dt.matFromBed(path2, structure2)
at.makeSymmetric(contacts1)
at.makeSymmetric(contacts2)

compartments1 = ca.get_compartments(contacts1, structure1, "binding_data/{}_{}_{}kb_active_coverage.bed".format(format_celltype(cell_type1), chrom, res_kb), True)
compartments2 = ca.get_compartments(contacts2, structure2, "binding_data/{}_{}_{}kb_active_coverage.bed".format(format_celltype(cell_type2), chrom, res_kb), True)	

gen_coords = structure1.getGenCoords()

dists = normalize(dists)
compartment_diffs = np.abs(compartments1 - compartments2)
compartment_diffs = normalize(compartment_diffs)

dist_peaks = sg.find_peaks_cwt(dists, np.arange(1, 20))

gen_coords = structure1.getGenCoords()

with open("{}_dist_peaks.bed".format(chrom), "w") as out:
	for peak in dist_peaks:
		out.write("\t".join((structure1.chrom.name, str(gen_coords[peak]), str(gen_coords[peak] + structure1.chrom.res), str(compartments1[peak]), str(compartments2[peak]))))
		out.write("\n")
	out.close()

diff_peaks = sg.find_peaks_cwt(compartment_diffs, np.arange(1, 20))

with open("{}_comp_peaks.bed".format(chrom), "w") as out:
	for peak in diff_peaks:
		out.write("\t".join((structure1.chrom.name, str(gen_coords[peak]), str(gen_coords[peak] + structure1.chrom.res))))
		out.write("\n")
	out.close()
