from statsmodels.stats.multitest import multipletests
import sys
import os

in_path = sys.argv[1]
prefix = in_path.split(".")[0]
res = int(sys.argv[2])

ps = []

with open(in_path) as in_file:
	for line in in_file:
		line = line.strip().split()
		if line[0] != "\"logFC\"":	#skip header
			ps.append(float(line[4]))
	in_file.close()

reject, qs, alphacSidak, alphacBonf = multipletests(ps, method="fdr_bh")

i = 0

out1 = open(prefix + "_loc1.bed", "w")
out2 = open(prefix + "_loc2.bed", "w")

with open(in_path) as in_file:
	for line in in_file:
		line = line.strip().split()
		if line[0] != "\"logFC\"":	
			loc_id = line[0].strip("\"").split(":")
			chrom = loc_id[0]
			loc1, loc2 = loc_id[1].split(",")
			if qs[i] < 0.01:
				out1.write("\t".join((chrom, loc1, str(int(loc1) + res), line[1])))
				out1.write("\n")
				out2.write("\t".join((chrom, loc2, str(int(loc2) + res))))
				out2.write("\n")
			i += 1
	in_file.close()

out1.close()
out2.close()

os.system("bedtools intersect -a %s_loc1.bed -b mappability.bed -wb > %s_loc1_mappability.bed"%(prefix, prefix))
os.system("bedtools intersect -a %s_loc2.bed -b mappability.bed -wb > %s_loc2_mappability.bed"%(prefix, prefix))
os.system("paste %s_loc1_mappability.bed %s_loc2_mappability.bed | awk '$8 > 0.8 && $15 > 0.8 {print $2\"\t\"$10\"\t\"$4}' > %s_sig.tsv"%(prefix, prefix, prefix))
