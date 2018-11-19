import numpy as np
from scipy import stats as st
import sys
from matplotlib import pyplot as plt

mat1 = np.loadtxt(sys.argv[1], dtype=object)
enrichments1 = np.array(mat1[:,6], dtype=float)
mat2 = np.loadtxt(sys.argv[2], dtype=object)
enrichments2 = np.array(mat2[:,6], dtype=float)
print st.ttest_ind(enrichments1, enrichments2)

x_int_size = 0.1
x_start = -x_int_size/5.
x_end = max((max(enrichments1), max(enrichments2)))

print max(enrichments2)

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
counts, bounds, patches = plt.hist(enrichments1)
y_int_size = 2000
y_start = y_int_size/5.
y_end = counts[0] - y_int_size/5.
plt.title("Background A compartment", fontsize=14)
plt.xlabel("Enhancer coverage", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=8)
plt.savefig("background_enhancer_coverage")
plt.show()

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
counts, bounds, patches = plt.hist(enrichments2)
y_int_size = 10
y_start = y_int_size/5.
y_end = counts[0] - y_int_size/5.
plt.title("Intra-A relocalized loci", fontsize=14)
plt.xlabel("Enhancer coverage", fontsize=14)
plt.axis([x_start, x_end, y_start, y_end], frameon=False)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=6)	
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=8)
plt.savefig("relocalization_enhancer_coverage")
plt.show()
