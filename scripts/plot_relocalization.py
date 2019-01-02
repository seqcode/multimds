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
res_kb = 32

max_dists = []
max_gencoords = []
min_gencoords = []

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
for strain in ("Scer", "Suva"):
	chrom_name = "{}_{}".format(strain, chrom_num)
	os.system("python ~/git/multimds/multimds.py --full -P 0.1 -w 0 ctrl_{}_{}kb.bed galactose_{}_{}kb.bed".format(chrom_name, res_kb, chrom_name, res_kb))
	struct1 = dt.structure_from_file("ctrl_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
	struct2 = dt.structure_from_file("galactose_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
	dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())]
	max_dists.append(max(dists))
	max_gencoords.append(max(struct1.getGenCoords()))
	min_gencoords.append(min(struct1.getGenCoords()))
	plt.plot(struct1.getGenCoords(), dists, label=strain, lw=4)

#define offsets
xmin = min(min_gencoords)
xmax = max(max_gencoords)
x_range = xmax - xmin
x_start = xmin - x_range/25.
x_end = xmax + x_range/25.

ymin = 0
ymax = max(max_dists)
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)

gen_coord = struct1.getGenCoords()[struct1.get_rel_index(gene_loc)]
plt.scatter([gen_coord], [0.005], c="g", s=50, marker="*")
plt.annotate(gene_name, (gen_coord+20000, 0.01))

plt.legend(frameon=False)
plt.savefig(gene_name)
plt.show()
