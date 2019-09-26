from matplotlib import pyplot as plt
import numpy as np
import os

gene_names = ("Hxt1", "Has1", "Tda1", "Gal1", "Gal7", "Gal10", "Gal3", "Gal4", "Gal2")
logfcs = np.zeros((len(gene_names), 1))

for i, gene_name in enumerate(gene_names):
	os.system("bedtools intersect -a ctrl_IP.bedgraph -b {}.bed > ctrl_{}_Nup60.bed".format(gene_name, gene_name))
	os.system("bedtools intersect -a galactose_IP.bedgraph -b {}.bed > galactose_{}_Nup60.bed".format(gene_name, gene_name))

	ctrl_counts = np.loadtxt("ctrl_{}_Nup60.bed".format(gene_name), usecols=range(1,4))
	weighted_ctrl_counts = (ctrl_counts[:,1] - ctrl_counts[:,0])*ctrl_counts[:,2]	#weight by length of interval
	galactose_counts = np.loadtxt("galactose_{}_Nup60.bed".format(gene_name), usecols=range(1,4))
	weighted_galactose_counts = (galactose_counts[:,1] - galactose_counts[:,0])*galactose_counts[:,2]

	logfcs[i] = np.log(np.mean(weighted_galactose_counts)/np.mean(weighted_ctrl_counts))

#no need to do anything fancy when defining our figure
fig, ax = plt.subplots()
plt.subplot2grid((10,10), (0,0), 10, 5, frameon=False)

plt.pcolor(logfcs, cmap=plt.cm.coolwarm, vmin=-1, vmax=1)

#plot ticks
indices = np.arange(len(logfcs)) + 0.5
labels = gene_names
plt.yticks(indices, labels)
plt.xticks(indices, [])
plt.tick_params(top=False, right=False, left=False, bottom=False, labelsize=12)		#don't want any ticks showing
cbaxes = fig.add_axes([0.3, 0.1, 0.02, 0.4]) 
plt.colorbar(cax=cbaxes)
plt.savefig("sup13a.svg")
