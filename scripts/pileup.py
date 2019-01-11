from matplotlib import pyplot as plt
import sys
import os
import numpy as np

gene = sys.argv[1]
os.system("bedtools intersect -a nup60_sig.bed -b {}.bed -f 0.5 > {}_Nup60_peak.bed".format(gene, gene))
os.system("bedtools intersect -a ctrl_IP.bedgraph -b {}.bed > ctrl_{}_Nup60.bed".format(gene, gene))
os.system("bedtools intersect -a galactose_IP.bedgraph -b {}.bed > galactose_{}_Nup60.bed".format(gene, gene))

plt.subplot2grid((10,10), (0,0), 9, 10, frameon=False)

max_num_tags = 0

first = True
with open ("ctrl_{}_Nup60.bed".format(gene)) as infile:
	for line in infile:
		line = line.strip().split()
		start = float(line[1])/1000
		end = float(line[2])/1000
		num_tags = int(line[3])
		if num_tags > max_num_tags:
			max_num_tags = num_tags
		if first:
			plt.plot([start,end], [num_tags,num_tags], c="b", label="glucose")
			first = False
		else:
			plt.plot([start,end], [num_tags,num_tags], c="b")
	infile.close()

first = True 

with open ("galactose_{}_Nup60.bed".format(gene)) as infile:
	for line in infile:
		line = line.strip().split()
		start = float(line[1])/1000
		end = float(line[2])/1000
		num_tags = int(line[3])
		if num_tags > max_num_tags:
			max_num_tags = num_tags
		if first:
			plt.plot([start,end], [num_tags,num_tags], c="g", label="galactose")
			first = False
		else:
			plt.plot([start,end], [num_tags,num_tags], c="g")
	infile.close()

peak = np.loadtxt("{}_Nup60_peak.bed".format(gene), dtype=object)
plt.plot([float(peak[1])/1000, float(peak[2])/1000], [30,30], c="k", label="Differential Nup60 peak")

plt.xlabel("Genomic coordinate (kb)", fontsize=12)
plt.ylabel("Tag count", fontsize=12)

#define offsets
gene_loc = np.loadtxt("{}.bed".format(gene), dtype=object)
gene_start = float(gene_loc[1])/1000
gene_end = float(gene_loc[2])/1000
xmin = gene_start
xmax = gene_end
x_range = xmax - xmin
x_start = xmin - x_range/25.
x_end = xmax + x_range/25.

ymin = 0
ymax = max_num_tags
y_range = ymax - ymin
y_start = ymin - y_range/25.
y_end = ymax + y_range/10.

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

plt.title(gene)
plt.legend(loc=2, fontsize=8, frameon=False, shadow=False)

plt.savefig("{}_Nup60".format(gene))
plt.show()