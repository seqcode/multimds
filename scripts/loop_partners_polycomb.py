import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as st
import sys

res_kb = int(sys.argv[1])

if os.path.isfile("polycomb_enrichment.txt"):
	os.system("rm polycomb_enrichment.txt")

if os.path.isfile("enhancer_enrichment.txt"):
	os.system("rm enhancer_enrichment.txt")

chroms = ["chr{}".format(chrom_num) for chrom_num in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)]

partners = {}
for chrom in chroms:
	partners[chrom] = {}

for chrom in chroms:
	with open("{}_{}kb_edgeR_output_sig.tsv".format(chrom, res_kb)) as infile:
		for line in infile:
			line = line.strip().split()
			loc1 = int(line[0])
			loc2 = int(line[1])
			fc = float(line[2])
			try:
				old_fc = partners[chrom][loc1][1]
				if np.abs(fc) > np.abs(old_fc):
					partners[chrom][loc1] = (loc2, fc)
			except KeyError:
				partners[chrom][loc1] = (loc2, fc)
			try:
				old_fc = partners[chrom][loc2][1]
				if np.abs(fc) > np.abs(old_fc):
					partners[chrom][loc2] = (loc1, fc)
			except KeyError:
				partners[chrom][loc2] = (loc1, fc)
		infile.close()

with open("peaks_filtered_GM12878_only_enhancer.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom = int(line[0])
		loc = line[1]
		try:
			partner, fc = partners[chrom][loc]
			if fc < 0:	#loop in K562 only
				os.system("cat binding_data/wgEncodeBroadHistoneK562H3k27me3StdPk_%dkb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> polycomb_enrichment.txt"%(res_kb, chrom, partner))
			else:	#loop in GM12878 only
				os.system("cat binding_data/GM12878_enhancers_%dkb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> enhancer_enrichment.txt"%(res_kb, chrom, partner))
		except KeyError:
			pass
	in_file.close()


with open("peaks_filtered_K562_only_enhancer.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom = int(line[0])
		loc = line[1]
		try:
			partner, fc = partners[chrom][loc]
			if fc > 0:	#loop in GM12878 only
				os.system("cat binding_data/wgEncodeBroadHistoneGm12878H3k27me3StdPkV2_%dkb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> polycomb_enrichment.txt"%(res_kb, chrom, partner))
			else:	#loop in K562 only
				os.system("cat binding_data/K562_enhancers_%dkb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> enhancer_enrichment.txt"%(res_kb, chrom, partner))
		except KeyError:
			pass
	in_file.close()

with open("peaks_filtered_both_enhancer.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom = int(line[0])
		loc = line[1]
		try:
			partner, fc = partners[chrom][loc]
			os.system("cat binding_data/GM12878_enhancers_%dkb_windows_enrichment.bed | awk '$1 == \"%s\" && $2 == %s {print $4}' >> polycomb_enrichment.txt"%(res_kb, chrom, partner))
		except KeyError:
			pass
	in_file.close()

os.system("bedtools coverage -a A_background_filtered.bed -b binding_data/wgEncodeBroadHistoneGm12878H3k27me3StdPkV2.broadPeak > A_background_filtered_polycomb.bed")

partner_enrichment = np.loadtxt("polycomb_enrichment.txt")
mat = np.loadtxt("A_background_filtered_polycomb.bed", dtype=object)
background_enrichment = np.array(mat[:,3], dtype=float)

print st.ttest_ind(background_enrichment, partner_enrichment)

plt.boxplot([background_enrichment, partner_enrichment], labels=("Background A compartment", "Loop partners"))
plt.ylabel("H3K27me3 enrichment")
plt.savefig("polycomb_enrichment")
