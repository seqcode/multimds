import sys
sys.path.append("..")
import data_tools as dt
import linear_algebra as la
from matplotlib import pyplot as plt
import numpy as np
import compartment_analysis as ca
from scipy import stats as st

cell_type1 = sys.argv[1]
cell_type2 = sys.argv[2]
res_kb = int(sys.argv[3])

struct1 = dt.structure_from_file("{}_21_{}kb_structure.tsv".format(cell_type1, res_kb))
struct2 = dt.structure_from_file("{}_21_{}kb_structure.tsv".format(cell_type2, res_kb))
gen_coords = np.array(struct1.getGenCoords())
dists = np.array([la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())])

mat1 = dt.matFromBed("hic_data/{}_21_{}kb.bed".format(cell_type1, res_kb), struct1)
comps1 = ca.get_compartments(mat1)
mat2 = dt.matFromBed("hic_data/{}_21_{}kb.bed".format(cell_type2, res_kb), struct2)
comps2 = ca.get_compartments(mat2)

r, p = st.pearsonr(comps1, comps2)
if r < 0:
	comps1 = -comps1

comp_diffs = np.abs(comps1 - comps2)

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.plot(gen_coords, dists/max(dists), lw=2, label="Relocalization", zorder=1)
plt.plot(gen_coords, comp_diffs/max(comp_diffs), lw=2, label="Compartment score change", zorder=1)

plt.title("{} vs {}".format(cell_type1, cell_type2), fontsize=15)

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
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.locator_params(axis="x", nbins=3)
plt.locator_params(axis="y", nbins=3)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=15)

gen_coord = 47400000
plt.scatter([gen_coord], [0.005], c="r", s=100, marker="o", zorder=2)

plt.legend(frameon=False)

plt.savefig("{}_{}_relocalization".format(cell_type1, cell_type2))
