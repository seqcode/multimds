import os
import sys
sys.path.append("/home/lur159/git/miniMDS")
import data_tools as dt
import linear_algebra as la
from matplotlib import pyplot as plt
import numpy as np

gene_name = sys.argv[1]
chrom_num = sys.argv[2]
gene_loc = int(sys.argv[3])
strain = sys.argv[4]
res_kb = 32

chrom_name = "{}_{}".format(strain, chrom_num)
os.system("python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_{}_{}kb.bed hic_data/galactose_{}_{}kb.bed".format(chrom_name, res_kb, chrom_name, res_kb))
struct1 = dt.structure_from_file("ctrl_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
struct2 = dt.structure_from_file("galactose_{}_{}kb_structure.tsv".format(chrom_name, res_kb))
gen_coords = np.array(struct1.getGenCoords())/1000
dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())]

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.plot(gen_coords, dists, lw=4)

plt.xlabel("chr{}".format(chrom_num), fontsize=20)

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
plt.locator_params(axis="x", nbins=3)
plt.locator_params(axis="y", nbins=3)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=15)

gen_coord = struct1.getGenCoords()[struct1.get_rel_index(gene_loc)]/1000
plt.scatter([gen_coord], [0.005], c="g", s=80, marker="*")
plt.annotate(gene_name, (gen_coord+20, 0.01), fontsize=16)

plt.savefig("{}_{}".format(strain, gene_name))
plt.show()
