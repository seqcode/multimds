cell_type1 = "GM12878_combined"
cell_type2 = "K562"
res = 100000

boundaries = []

with open("{}_tadlib_output.txt".format(cell_type1)) as in_file:
	for line in in_file:
		line = line.split()
		boundary1 = line[0] + "-" + line[1]
		if boundary1 not in boundaries:
			boundaries.append(boundary1)
		boundary2 = line[0] + "-" + line[2]
		if boundary2 not in boundaries:
			boundaries.append(boundary2)
	in_file.close()

unique = []

with open("{}_tadlib_output.txt".format(cell_type2)) as in_file:
	for line in in_file:
		line = line.split()
		boundary1 = line[0] + "-" + line[1]
		if boundary1 not in boundaries and boundary1 not in unique:
			unique.append(boundary1)
		boundary2 = line[0] + "-" + line[2]
		if boundary2 not in boundaries and boundary2 not in unique:
			unique.append(boundary2)
	in_file.close()

with open("{}_{}_{}kb_differential_tad_boundaries.bed".format(cell_type1, cell_type2, res/1000), "w") as out_file:
	for boundary in unique:
		chrom, loc = boundary.split("-")
		out_file.write("\t".join(("chr{}".format(chrom), loc, str(int(loc) + res))))
		out_file.write("\n")
	out_file.close()
