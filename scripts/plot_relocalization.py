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
prefix1 = sys.argv[4]
prefix2 = sys.argv[5]
res_kb = 32

max_dists = []
max_gencoords = []

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
for strain in ("Scer", "Suva"):
	chrom_name = "{}_{}".format(strain, chrom_num)
	os.system("python ~/git/multimds/multimds.py --full -P 0.1 -w 0 {}_{}_{}kb.bed {}_{}_{}kb.bed".format(prefix1, chrom_name, res_kb, prefix2, chrom_name, res_kb))
	struct1 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(prefix1, chrom_name, res_kb))
	struct2 = dt.structure_from_file("{}_{}_{}kb_structure.tsv".format(prefix2, chrom_name, res_kb))
	dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())]
	max_dists.append(max(dists))
	max_gencoords.append(max(struct1.getGenCoords()))
	plt.plot(struct1.getGenCoords(), dists, label=strain, lw=4)

x_int_size = 200000
ys = dists
y_int_size = 0.01
x_start = -x_int_size/4.
x_end = max(max_gencoords) + x_int_size/5.
y_start = -y_int_size/5.
y_end = max(max_dists) + y_int_size/5.

plt.title("chr{}".format(chrom_num), fontsize=14)
plt.xlabel("Genomic coordinate", fontsize=14)
plt.ylabel("Relocalization", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end],frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)	 
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)
gen_coord = struct1.getGenCoords()[struct1.get_rel_index(gene_loc)]
plt.scatter([gen_coord], [0.005], c="g", s=50, marker="*")
plt.annotate(gene_name, (gen_coord+20000, 0.005))
plt.legend()
plt.show()
#plt.savefig(gene_name)
