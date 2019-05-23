import sys
sys.path.append("..")
import compartment_analysis as ca
import data_tools as dt
from scipy import stats as st
from matplotlib import pyplot as plt
import numpy as np
import linear_algebra as la
from scipy import signal as sg
import os

chrom = 19
res_kb = 100

os.system("python ../multimds.py -o test_ hic_data/GM12878_combined_{}_{}kb.bed hic_data/K562_{}_{}kb.bed".format(chrom, res_kb, chrom, res_kb))

struct1 = dt.structure_from_file("test_GM12878_combined_{}_{}kb_structure.tsv".format(chrom, res_kb))	
struct2 = dt.structure_from_file("test_K562_{}_{}kb_structure.tsv".format(chrom, res_kb))

mat1 = dt.matFromBed("hic_data/GM12878_combined_{}_{}kb.bed".format(chrom, res_kb), struct1)
comps1 = ca.get_compartments(mat1)
mat2 = dt.matFromBed("hic_data/K562_{}_{}kb.bed".format(chrom, res_kb), struct2)
comps2 = ca.get_compartments(mat2)

r, p = st.pearsonr(comps1, comps2)
if r < 0:
	comps1 = -comps1

comp_diffs = np.abs(comps1 - comps2)

dists = np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())])
dist_peaks = sg.find_peaks_cwt(dists, np.arange(1,10))

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
gen_coords = struct1.getGenCoords()
plt.plot(gen_coords, comp_diffs/max(comp_diffs), lw=2, color=(0.75,0,0), label="Compartment score change", zorder=1)
plt.plot(gen_coords, dists/max(dists), lw=2, color=(0,0,0.75), label="Relocalization", zorder=1)

for dist_peak in dist_peaks:
	if comp_diffs[dist_peak] < 0.2:
		plt.scatter([gen_coords[dist_peak]], [0.1], color=(0,0,1), s=40, zorder=2)
	else:
		plt.scatter([gen_coords[dist_peak]], [0.1], color=(0.25,0.25,0.25), s=40, zorder=2)


plt.xlabel("Genomic coordinate", fontsize=20)
plt.ylabel("Normalized change", fontsize=20)

#define offsets
xmin = min(gen_coords)
xmax = max(gen_coords)
x_range = xmax - xmin
x_start = xmin - x_range/25.
x_end = xmax + x_range/25.

ymin = 0
ymax = 1
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=2)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=18)

plt.plot([x_start, x_end], [0.2, 0.2], color=(0.5,0.5,0.5), linestyle="--")

plt.legend(frameon=False, loc=2, fontsize=16)

plt.savefig("sup16")
