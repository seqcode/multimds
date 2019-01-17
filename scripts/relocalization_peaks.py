import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import compartment_analysis as ca
import os
import linear_algebra as la
from scipy import signal as sg
from hmmlearn import hmm

def normalize(values):
	return np.array(values)/max(values)

def format_celltype(cell_type):
	if cell_type == "KBM7":
		return "K562"	#substitute
	else:
		formatted = cell_type.split("_")[0]
		return formatted[0].upper() + formatted[1:len(formatted)].lower()

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

cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
chrom = sys.argv[3]
#centromere = sys.argv[4]
#num_partitions = sys.argv[5]
smoothing_parameter = float(sys.argv[6])
res = int(sys.argv[7])
res_kb = res/1000
#n = 1

#path1 = "hic_data/{}_{}_{}kb_filtered.bed".format(cell_type1, chrom, res_kb)
#path2 = "hic_data/{}_{}_{}kb_filtered.bed".format(cell_type2, chrom, res_kb)
path1 = "hic_data/{}_{}_{}kb.bed".format(cell_type1, chrom, res_kb)
path2 = "hic_data/{}_{}_{}kb.bed".format(cell_type2, chrom, res_kb)

#min_error = sys.float_info.max
#for iteration in range(n):
	#os.system("python ../multimds.py -m {} -N {} -o {}_ {} {}".format(centromere, num_partitions, iteration, path1, path2))

#os.system("python /home/lur159/git/multimds/multimds.py --full {} {}".format(path1, path2))
		
#load structures
#structure1 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
#structure2 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type2, chrom, res_kb))
structure1 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type1, chrom, res_kb))	
structure2 = dt.structure_from_file("hic_data/{}_{}_{}kb_structure.tsv".format(cell_type2, chrom, res_kb))

#rescale
structure1.rescale()
structure2.rescale()

#make structures compatible
dt.make_compatible((structure1, structure2))

#align
#r, t = la.getTransformation(structure1, structure2)
#structure1.transform(r,t)

	#calculate error
	#coords1 = np.array(structure1.getCoords())
	#coords2 = np.array(structure2.getCoords())
	#error = np.mean([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)])
	#if error < min_error:
	#	min_error = error
	#	best_iteration = iteration

#for iteration in range(n):
#	if iteration == best_iteration:
		#load structures
#		structure1 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
#		structure2 = dt.structure_from_file("/data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type2, chrom, res_kb))
#	else:
#		os.system("rm /data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type1, chrom, res_kb))	
#		os.system("rm /data/drive1/test/archive/multimds/scripts/hic_data/{}_{}_{}_{}kb_filtered_structure.tsv".format(iteration, cell_type2, chrom, res_kb))	

#rescale
#structure1.rescale()
#structure2.rescale()

#make structures compatible
#dt.make_compatible((structure1, structure2))

#align
#r, t = la.getTransformation(structure1, structure2)
#structure1.transform(r,t)

#calculate changes
coords1 = np.array(structure1.getCoords())
coords2 = np.array(structure2.getCoords())
dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(coords1, coords2)]

#compartments
chrom1 = dt.chromFromBed(path1)
chrom2 = dt.chromFromBed(path2)
chrom1.res = 100000	#reduce res to reduce RAM usage in compartment calculation
chrom2.res = 100000
chrom1.minPos = int(np.floor(float(chrom1.minPos)/chrom1.res)) * chrom1.res	#round
chrom1.maxPos = int(np.ceil(float(chrom1.maxPos)/chrom1.res)) * chrom1.res
chrom2.minPos = int(np.floor(float(chrom2.minPos)/chrom2.res)) * chrom2.res	#round
chrom2.maxPos = int(np.ceil(float(chrom2.maxPos)/chrom2.res)) * chrom2.res

low_struct1 = dt.structureFromBed(path1, chrom1)
low_struct2 = dt.structureFromBed(path2, chrom2)
dt.make_compatible((low_struct1, low_struct2))
contacts1 = dt.matFromBed(path1, low_struct1)		
contacts2 = dt.matFromBed(path2, low_struct2)

enrichments = np.loadtxt("binding_data/{}_{}_100kb_active_coverage.bed".format(format_celltype(cell_type1), chrom), usecols=6)
bin_nums = low_struct1.nonzero_abs_indices() + low_struct1.chrom.minPos/low_struct1.chrom.res
enrichments = enrichments[bin_nums]
compartments1 = np.array(ca.get_compartments(contacts1, enrichments))

enrichments = np.loadtxt("binding_data/{}_{}_100kb_active_coverage.bed".format(format_celltype(cell_type2), chrom), usecols=6)
bin_nums = low_struct2.nonzero_abs_indices() + low_struct2.chrom.minPos/low_struct2.chrom.res
enrichments = enrichments[bin_nums]
compartments2 = np.array(ca.get_compartments(contacts2, enrichments))

gen_coords = structure1.getGenCoords()

#dists = normalize(dists)
compartment_diffs = np.abs(compartments1 - compartments2)
#compartment_diffs = normalize(compartment_diffs)

#smoothed_dists = sg.cwt(dists, sg.ricker, [smoothing_parameter])[0]
#dist_peaks = call_peaks(smoothed_dists)
dist_peaks = sg.find_peaks_cwt(dists, np.arange(1,10))
#randomize
#num_peaks = len(dist_peaks)
#dist_peaks = np.random.randint(len(dists), size=num_peaks)

#smoothed_diffs = sg.cwt(compartment_diffs, sg.ricker, [smoothing_parameter])[0]
#diff_peaks = call_peaks(smoothed_diffs)
#diff_peaks = sg.find_peaks_cwt(compartment_diffs, np.arange(1,10))
#randomize
#num_peaks = len(diff_peaks)
#diff_peaks = np.random.randint(len(compartment_diffs), size=num_peaks)

high_coords = structure1.getGenCoords()
low_coords = low_struct1.getGenCoords()

high_comps1 = []
high_comps2 = []

for high_coord in high_coords:
	low_index1 = low_struct1.get_rel_index(high_coord)
	assert low_index1 is not None
	low_comp1 = compartments1[low_index1]
	high_comps1.append(low_comp1)
	low_index2 = low_struct2.get_rel_index(high_coord)
	assert low_index2 is not None
	low_comp2 = compartments2[low_index2]
	high_comps2.append(low_comp2)

with open("{}_dist_peaks.bed".format(chrom), "w") as out:
	for peak in dist_peaks:
		#start, end = peak
		#peak_dists = dists[start:end]
		#max_dist_index = np.argmax(peak_dists) + start
		max_dist_index = peak
		#out.write("\t".join(("{}".format(structure1.chrom.name), str(gen_coords[start]), str(gen_coords[end]), str(gen_coords[max_dist_index]))))
		out.write("\t".join((structure1.chrom.name, str(high_coords[max_dist_index]), str(high_coords[max_dist_index] + structure1.chrom.res), str(high_comps1[max_dist_index]), str(high_comps2[max_dist_index]))))
		#out.write("\t".join((structure1.chrom.name, str(high_coords[start]), str(high_coords[end]), str(high_comps1[max_dist_index]), str(high_comps2[max_dist_index]))))
		out.write("\n")
	out.close()
sys.exit(0)

with open("{}_comp_peaks.bed".format(chrom), "w") as out:
	for peak in diff_peaks:
		#start, end = peak
		#peak_diffs = compartment_diffs[start:end]
		#max_diff_index = np.argmax(peak_diffs) + start
		#max_diff_index = peak
		#out.write("\t".join((structure1.chrom.name, str(high_coords[max_diff_index]), str(high_coords[max_diff_index] + low_struct1.chrom.res))))
		out.write("\t".join((structure1.chrom.name, str(gen_coords[peak]), str(gen_coords[peak] + structure1.chrom.res))))
		#out.write("\t".join((structure1.chrom.name, str(gen_coords[start]), str(gen_coords[end]))))
		out.write("\n")
	out.close()
