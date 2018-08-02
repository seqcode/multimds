import numpy as np
from scipy import stats as st
import sys
from matplotlib import pyplot as plt

mat1 = np.loadtxt(sys.argv[1], dtype=object)
enrichments1 = np.array(mat1[:,6], dtype=float)
mat2 = np.loadtxt(sys.argv[2], dtype=object)
enrichments2 = np.array(mat2[:,6], dtype=float)
print st.ttest_ind(enrichments1, enrichments2)
plt.hist(enrichments1)
plt.savefig("relocalization_H3K27ac_coverage")
plt.show()
plt.hist(enrichments2)
plt.savefig("background_H3K27ac_coverage")
plt.show()
