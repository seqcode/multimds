import sys
sys.path.append("..")
import data_tools as dt
from matplotlib import pyplot as plt
import numpy as np

res = 100000
res_kb = res/1000
loc = 47400000

for celltype in ("K562", "HMEC", "HUVEC", "IMR90"):
	path = "hic_data/{}_21_{}kb.bed".format(celltype, res_kb)
	struct = dt.structureFromBed(path)
	mat = dt.matFromBed(path, struct)
	gen_coords = struct.getGenCoords()
	filtered_gen_coords = []
	filtered_counts = []
	counts = mat[struct.get_rel_index(loc)]
	start = struct.get_rel_index(43000000)
	end = struct.get_rel_index(47300000)
	counts = counts[start:end]
	gen_coords = gen_coords[start:end]
	plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
	plt.plot(gen_coords, counts, lw=4)

	plt.title(celltype, fontsize=12)

	#define offsets
	xmin = min(gen_coords)
	xmax = max(gen_coords)
	x_range = xmax - xmin
	x_start = xmin - x_range/25.
	x_end = xmax + x_range/25.

	ymin = 0
	ymax = max(counts)
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
	plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=10)

	gen_coord = 46900000
	plt.plot([gen_coord, gen_coord], [ymin, ymax], c="g", linestyle="--", zorder=2)

	plt.savefig("{}_{}kb_4C".format(celltype, res_kb))
