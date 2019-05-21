from matplotlib import pyplot as plt
import seaborn as sns
import os
import pandas as pd

os.system("paste peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed peaks_filtered_{}_coverage.bed | cut -f 7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126 > peaks_filtered_all_coverage.bed".format("GM12878_H3K27ac", "GM12878_H3K4me1", "GM12878_H3K4me3", "GM12878_H3K9ac", "GM12878_H2AZ", "GM12878_H3K4me2", "GM12878_H3K27me3", "K562_H3K27ac", "K562_H3K4me1", "K562_H3K4me3", "K562_H3K9ac", "K562_H2AZ", "K562_H3K4me2", "K562_H3K27me3"))

sns.set(font_scale=0.5)
all_coverage = pd.read_table("peaks_filtered_all_coverage.bed", names=("GM12878 H3K27ac", "GM12878 H3K4me1", "GM12878 H3K4me3", "GM12878 H3K9ac", "GM12878 H2AZ", "GM12878 H3K4me2", "GM12878 H3K27me3", "K562 H3K27ac", "K562 H3K4me1", "K562 H3K4me3", "K562 H3K9ac", "K562 H2AZ", "K562 H3K4me2", "K562 H3K27me3"))
sns.clustermap(all_coverage, cmap="Reds", col_cluster=False, cbar_kws={"label": "Coverage"})
plt.savefig("hierarchical_clustering")
plt.show()
