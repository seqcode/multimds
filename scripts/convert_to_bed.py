import os

chrom_bins = {}

with open("GSE88952_Sc_Su.32000.bed") as in_file:
	for line in in_file:
		line = line.strip().split()
		chrom_bins[line[3]] = "{}\t{}\t{}".format(line[0], line[1], line[2])
	in_file.close()

if not os.path.isfile("ctrl_32kb.bed"):
	with open("ctrl_32kb.bed", "w") as out_file:
		with open("ctrl_32kb_matrix.txt") as in_file:
			for line in in_file:
				line = line.strip().split()
				bin1 = line[0]
				chrom_string1 = chrom_bins[bin1]
				bin2 = line[1]
				chrom_string2 = chrom_bins[bin2]
				if float(line[3]) != 0:
					out_file.write("\t".join((chrom_string1, chrom_string2, line[3])))
					out_file.write("\n")
			in_file.close()
		out_file.close()

if not os.path.isfile("galactose_32kb.bed"):
	with open("galactose_32kb.bed", "w") as out_file:
		with open("galactose_32kb_matrix.txt") as in_file:
			for line in in_file:
				line = line.strip().split()
				bin1 = line[0]
				chrom_string1 = chrom_bins[bin1]
				bin2 = line[1]
				chrom_string2 = chrom_bins[bin2]
				if float(line[3]) != 0:
					out_file.write("\t".join((chrom_string1, chrom_string2, line[3])))
					out_file.write("\n")
			in_file.close()
		out_file.close()
