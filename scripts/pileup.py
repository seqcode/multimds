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
			plt.plot([start,end], [num_tags,num_tags], c="c", label="glucose")
			first = False
		else:
			plt.plot([start,end], [num_tags,num_tags], c="c")
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
y_end = ymax + y_range/5.

#define axes with offsets
plt.axis([x_start, x_end, y_start, y_end], frameon=False)

#plot axes (black with line width of 4)
plt.axvline(x=x_start, color="k", lw=4)
plt.axhline(y=y_start, color="k", lw=4)

#plot ticks
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=5, labelsize=12)

peaks = np.loadtxt("{}_Nup60_peak.bed".format(gene), usecols=(1,2))
if len(peaks) > 0:
	if len(peaks.shape) == 1:
		peak = peaks
		plt.plot([peak[0]/1000, peak[1]/1000], [y_range/70., y_range/70.], c="r", label="Differential Nup60 peak")
	else:
		num_peaks = peaks.shape[1]
		for i in range(num_peaks):
			peak = peaks[i]
			if i == 0:
				plt.plot([peak[0]/1000, peak[1]/1000], [y_range/70., y_range/70.], c="r", label="Differential Nup60 peak")
			else:
				plt.plot([peak[0]/1000, peak[1]/1000], [y_range/70., y_range/70.], c="r")

with open ("{}_transcription_direction.bed".format(gene)) as infile:
	for line in infile:
		line = line.strip().split()
		start = float(line[1])/1000
		end = float(line[2])/1000
		plt.arrow(start, y_start + y_range/10., end-start, 0, facecolor="k", head_width=6, head_length=0.1)
		plt.annotate(line[3], (np.mean((start,end)), y_start + y_range/8.), fontsize=16)
	infile.close()

plt.legend(loc=2, fontsize=8, frameon=False, shadow=False)

plt.savefig("{}_Nup60".format(gene))
plt.show()
