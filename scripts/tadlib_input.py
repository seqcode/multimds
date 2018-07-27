import sys
sys.path.append("..")
import data_tools as dt
import os

cell_type = sys.argv[1]

os.system("mkdir -p {}_tadlib_input".format(cell_type))

for chrom in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22):
	path = "hic_data/{}_{}_100kb.bed".format(cell_type, chrom)
	structure = dt.structureFromBed(path)
	mat = dt.matFromBed(path, structure)
	points = structure.getPoints()
	with open("{}_tadlib_input/chr{}.txt".format(cell_type, chrom), "w") as out:
		for i in range(len(mat)):
			point_num1 = points[i].num
			for j in range(i):
				if mat[i,j] != 0:
					point_num2 = points[j].num
					out.write("\t".join((str(point_num1), str(point_num2), str(mat[i,j]))))
					out.write("\n")
		out.close()
