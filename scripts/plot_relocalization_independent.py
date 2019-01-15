import os
import sys
sys.path.append("/home/lur159/git/miniMDS")
import data_tools as dt
import linear_algebra as la
from matplotlib import pyplot as plt
import numpy as np

chrom_num = sys.argv[1]
gene_name = sys.argv[2]
gene_loc = int(sys.argv[3])
strain = sys.argv[4]
res_kb = 32

chrom_name = "{}_{}".format(strain, chrom_num)
os.system("python ../minimds.py -w 0 -o independent_ctrl_{}_{}kb_structure.tsv ctrl_{}_{}kb.bed".format(chrom_name, res_kb, chrom_name, res_kb))
os.system("python ../minimds.py -w 0 -o independent_galactose_{}_{}kb_structure.tsv galactose_{}_{}kb.bed".format(chrom_name, res_kb, chrom_name, res_kb))
struct1 = dt.structure_from_file("independent_ctrl_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
struct2 = dt.structure_from_file("independent_galactose_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
dt.make_compatible((struct1, struct2))
struct1.rescale()
struct2.rescale()
r, t = la.getTransformation(struct1, struct2)
struct1.transform(r,t)
gen_coords = struct1.getGenCoords()
dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())]
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.plot(gen_coords, dists, lw=4)

plt.xlabel("Genomic coordinate", fontsize=11)
plt.ylabel("Glucose vs. galactose relocalization", fontsize=11)

#define offsets
xmin = min(gen_coords)
xmax = max(gen_coords)
x_range = xmax - xmin
x_start = xmin - x_range/25.
x_end = xmax + x_range/25.

ymin = 0
ymax = max(dists)
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=10)

gen_coord = struct1.getGenCoords()[struct1.get_rel_index(gene_loc)]
plt.scatter([gen_coord], [0.005], c="g", s=50, marker="*")
plt.annotate(gene_name, (gen_coord+20000, 0.01))

plt.savefig("{}_{}_independent".format(strain, gene_name))
plt.show()
