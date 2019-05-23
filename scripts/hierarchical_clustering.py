from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(font_scale=0.5)
all_coverage = pd.read_table("peaks_filtered_all_coverage.bed", names=("GM12878 H3K27ac", "GM12878 H3K4me1", "GM12878 H3K4me3", "GM12878 H3K9ac", "GM12878 H2AZ", "GM12878 H3K4me2", "GM12878 H3K27me3", "K562 H3K27ac", "K562 H3K4me1", "K562 H3K4me3", "K562 H3K9ac", "K562 H2AZ", "K562 H3K4me2", "K562 H3K27me3"))
sns.clustermap(all_coverage, cmap="Reds", col_cluster=False, cbar_kws={"label": "Coverage"})
plt.savefig("fig10")
