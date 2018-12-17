import numpy as np
from scipy import stats as st
import sys
from matplotlib import pyplot as plt

mat1 = np.loadtxt(sys.argv[1], dtype=object)
enrichments1 = np.array(mat1[:,6], dtype=float)
mat2 = np.loadtxt(sys.argv[2], dtype=object)
enrichments2 = np.array(mat2[:,6], dtype=float)
print st.ttest_ind(enrichments1, enrichments2)

xs = enrichments1 
#need to know bins to get y range
bins = plt.hist(xs)	
plt.close()

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.xlabel("GM12878 enhancer coverage", fontsize=14)
plt.title("Relocalized", fontsize=14)

#define offsets
xmin = min(xs)
xmax = max(xs)
x_range = xmax - xmin
x_start = xmin - x_range/25.	#bigger offset for bar plot
x_end = xmax + x_range/25.

ymin = 0
ymax = max(bins[0])
y_range = ymax - ymin
#y_start = ymin - y_range/25.
y_start = 0
y_end = ymax + y_range/25.

#plot
plt.hist(xs, rwidth=0.8, bottom=y_start)

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.savefig("relocalization_enhancer_coverage")
plt.close()

xs = enrichments2 
#need to know bins to get y range
bins = plt.hist(xs)	
plt.close()

#start with a frameless plot (extra room on the left)
plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

#label axes
plt.xlabel("GM12878 enhancer coverage", fontsize=14)
plt.title("Background", fontsize=14)

#define offsets
xmin = min(xs)
xmax = max(xs)
x_range = xmax - xmin
x_start = xmin - x_range/25.	#bigger offset for bar plot
x_end = xmax + x_range/25.

ymin = 0
ymax = max(bins[0])
y_range = ymax - ymin
#y_start = ymin - y_range/25.
y_start = 0
y_end = ymax + y_range/25.

#plot
plt.hist(xs, rwidth=0.8, bottom=y_start)

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.savefig("background_enhancer_coverage")
plt.close()
