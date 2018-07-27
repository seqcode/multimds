import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as st

os.system("rm polycomb_enrichment.txt enhancer_enrichment.txt")

with open("peaks_filtered_GM12878_only_enhancer.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom = line[0]
		loc = line[1]
		os.system("cat %s_edgeR_output_sig.tsv | awk '$1 == %s || $2 == %s {print $0}' > partners.tsv"%(chrom, loc, loc))
		mat = np.loadtxt("partners.tsv", dtype=object)
		if len(mat) > 0:
			try:
				fcs = np.abs(np.array(mat[:,2], dtype=float))
				best_line = mat[np.argmax(fcs)]
			except IndexError:
				best_line = mat
			if best_line[0] == loc:
				partner = best_line[1]
			else:
				partner = best_line[0]
			fc = float(best_line[2])
			if fc < 0:	#loop in K562 only
				os.system("cat binding_data/wgEncodeBroadHistoneK562H3k27me3StdPk_100kb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> polycomb_enrichment.txt"%(chrom, partner))
			else:	#loop in GM12878 only
				os.system("cat GM12878_enhancers_100kb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> enhancer_enrichment.txt"%(chrom, partner))
	in_file.close()


with open("peaks_filtered_K562_only_enhancer.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom = line[0]
		loc = line[1]
		os.system("cat %s_edgeR_output_sig.tsv | awk '$1 == %s || $2 == %s {print $0}' > partners.tsv"%(chrom, loc, loc))
		mat = np.loadtxt("partners.tsv", dtype=object)
		if len(mat) > 0:
			try:
				fcs = np.abs(np.array(mat[:,2], dtype=float))
				best_line = mat[np.argmax(fcs)]
			except IndexError:
				best_line = mat
			if best_line[0] == loc:
				partner = best_line[1]
			else:
				partner = best_line[0]
			os.system("cat mappability.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' > mappability.txt"%(chrom, partner))
			mappability = np.loadtxt("mappability.txt")
			fc = float(best_line[2])
			if fc > 0:	#loop in GM12878 only
				os.system("cat binding_data/wgEncodeBroadHistoneGm12878H3k27me3StdPkV2_100kb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> polycomb_enrichment.txt"%(chrom, partner))
			else:	#loop in K562 only
				os.system("cat K562_enhancers_100kb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> enhancer_enrichment.txt"%(chrom, partner))
	in_file.close()

with open("peaks_filtered_both_enhancer.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom = line[0]
		loc = line[1]
		os.system("cat %s_edgeR_output_sig.tsv | awk '$2 == %s || $3 == %s {print $0}' > partners.tsv"%(chrom, loc, loc))
		mat = np.loadtxt("partners.tsv", dtype=object)
		if len(mat) > 0:
			try:
				fcs = np.abs(np.array(mat[:,2], dtype=float))
				best_line = mat[np.argmax(fcs)]
			except IndexError:
				best_line = mat
			if best_line[0] == loc:
				partner = best_line[1]
			else:
				partner = best_line[0]
			os.system("cat mappability.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' > mappability.txt"%(chrom, partner))
			mappability = np.loadtxt("mappability.txt")
			os.system("cat binding_data/GM12878_enhancers_100kb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> polycomb_enrichment.txt"%(chrom, partner))
	in_file.close()

os.system("bedtools coverage -a A_background_filtered.bed -b binding_data/wgEncodeBroadHistoneGm12878H3k27me3StdPkV2.broadPeak > A_background_filtered_polycomb.bed")

partner_enrichment = np.loadtxt("polycomb_enrichment.txt")
mat = np.loadtxt("A_background_filtered_polycomb.bed", dtype=object)
background_enrichment = np.array(mat[:,3], dtype=float)

print st.ttest_ind(background_enrichment, partner_enrichment)

plt.boxplot([background_enrichment, partner_enrichment], labels=("Background A compartment", "Loop partners"))
plt.ylabel("H3K27me3 enrichment")
plt.savefig("polycomb_enrichment")
