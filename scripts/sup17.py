import numpy as np
from matplotlib import pyplot as plt

peak_comps = np.loadtxt("peaks_filtered.bed", usecols=(3,4))
background_comps = np.loadtxt("A_background_filtered.bed", usecols=(3,4))

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("Compartment change", fontsize=12)

#define offsets
ys = [background_comps[:,0] - background_comps[:,1], peak_comps[:,0] - peak_comps[:,1], ]
xs = range(len(ys))

xmin = 0
xmax = 0.4
x_range = xmax - xmin
x_start = xmin - x_range/2.	#larger offset for boxplot
x_end = xmax + x_range/10.

ymin = min([min(y) for y in ys])
ymax = max([max(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, positions=(0, 0.3), notch=True, patch_artist=True, labels=("Background", "Relocalized"), medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=1, labelsize=9)

plt.savefig("sup17a")
plt.close()

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.ylabel("Compartment score", fontsize=12)

#define offsets
ys = [background_comps[:,0], background_comps[:,1], peak_comps[:,0], peak_comps[:,1], ]
xs = range(len(ys))

xmin = 0
xmax = 1.5
x_range = xmax - xmin
x_start = xmin - x_range/5.	#larger offset for boxplot
x_end = xmax + x_range/10.

ymin = min([min(y) for y in ys])
ymax = max([max(y) for y in ys])
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/25.

#plot data
plt.boxplot(ys, positions=(0, 0.5, 1, 1.5), notch=True, patch_artist=True, labels=("GM12878\nbackground", "K562\nbackground", "GM12878\nrelocalized", "K562\nrelocalized"), medianprops=dict(linestyle="none"))	#boxplot has built-in support for labels, unlike barplot

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=10)

plt.savefig("sup17b")
