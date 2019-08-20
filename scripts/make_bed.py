import sys

num = sys.argv[1]

i = 0	#index
with open("sim{}_chr21.bed".format(num), "w") as out_file:
	with open("sim{}_chr21.fastq".format(num)) as in_file:
		for line in in_file:
			if line[0] == "@":	#header lines only
				line = line.strip().split()
				if i%2 == 1 and line[2] == "HIC":
					chrom1, start1 = line[3].split(":")
					chrom2, start2 = line[4].split(":")
					out_file.write("\t".join((chrom1, start1, str(int(start1)+1), chrom2, start2, str(int(start2)+1), "1")))
					out_file.write("\n")
				i += 1
		in_file.close()
	out_file.close()
