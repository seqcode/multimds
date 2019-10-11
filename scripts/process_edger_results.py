with open("nup60_sig.bed", "w") as outfile:	
	with open("nup60_edgeR_results.tsv") as infile:	
		for line in infile:	
			line = line.strip().split()	
			if line[0] != "\"logFC\"":	
				if float(line[4]) < 0.01:	
					chrom, bounds = line[0].strip("\"").split(":")	
					start, end = bounds.split("-")	
					outfile.write("\t".join((chrom, start, end)))	
					outfile.write("\n")	
		infile.close()	
	outfile.close()
