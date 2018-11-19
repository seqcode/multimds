import numpy as np
from scipy import stats as st
import sys
from matplotlib import pyplot as plt

mat1 = np.loadtxt(sys.argv[1], dtype=object)
enrichments1 = np.array(mat1[:,6], dtype=float)
mat2 = np.loadtxt(sys.argv[2], dtype=object)
enrichments2 = np.array(mat2[:,6], dtype=float)
print st.ttest_ind(enrichments1, enrichments2)

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.title("Intra-A relocalized loci")
plt.xlabel("Enhancer coverage")
bins = plt.hist(enrichments1)
print bins
#plt.axis([x_start, x_end, y_start, y_end], frameon=False)
#plt.axvline(x=x_start, color="k", lw=4)
#plt.axhline(y=y_start, color="k", lw=6)
plt.savefig("relocalization_H3K27ac_coverage")
plt.show()

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)
plt.title("Background A compartment")
plt.xlabel("Enhancer coverage")
bins = plt.hist(enrichments2)
print bins
plt.savefig("background_H3K27ac_coverage")
plt.show()
