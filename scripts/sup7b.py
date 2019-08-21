import matplotlib
matplotlib.use("Agg")
from multimds import data_tools as dt
from matplotlib import pyplot as plt
from multimds import linear_algebra as la
import os

os.system("python minimds.py sim1_chr21_100kb.bed")
os.system("python minimds.py sim2_chr21_100kb.bed")
struct1 = dt.structure_from_file("sim1_chr21_100kb_structure.tsv")
struct2 = dt.structure_from_file("sim2_chr21_100kb_structure.tsv")
struct1.rescale()
struct2.rescale()
r,t = la.getTransformation(struct1, struct2)
struct1.transform(r,t)

gen_coords = struct1.getGenCoords()
dists = [la.calcDistance(coord1, coord2) for coord1, coord2 in zip(struct1.getCoords(), struct2.getCoords())]

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.plot(gen_coords, dists, lw=2)

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

plt.xlabel("Genomic coordinate", fontsize=12)
plt.ylabel("Relocalization distance", fontsize=12)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=10)

#plt.show()
plt.savefig("sup7a")
