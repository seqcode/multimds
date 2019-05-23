from matplotlib import pyplot as plt
import numpy as np

gene_names = ("Hxt1", "Has1", "Tda1", "Gal1", "Gal7", "Gal10", "Gal3", "Gal4", "Gal2")
sgds = np.array(["YHR094C", "YMR290C", "YMR291W", "YBR020W", "YBR018C", "YBR019C", "YDR009W", "YPL248C", "YLR081W"])
logfcs = np.zeros((len(sgds), 1))

with open("rnaseq_counts.tsv") as infile:
	for line in infile:
		line = line.strip().split()
		if line[0] in sgds:
			index = np.where(sgds == line[0])[0][0]
			fc = np.mean([float(line[i]) for i in range(4,7)])/np.mean([float(line[i]) for i in range(1,4)])
			logfcs[index][0] = np.log(fc)
	infile.close()

#no need to do anything fancy when defining our figure
fig, ax = plt.subplots()
plt.subplot2grid((10,10), (0,0), 10, 5, frameon=False)

plt.pcolor(logfcs, cmap=plt.cm.coolwarm, vmin=-8, vmax=8)

#plot ticks
indices = np.arange(len(logfcs)) + 0.5
labels = gene_names
plt.yticks(indices, labels)
plt.xticks(indices, [])
plt.tick_params(top=False, right=False, left=False, bottom=False, labelsize=12)		#don't want any ticks showing
cbaxes = fig.add_axes([0.3, 0.1, 0.02, 0.4]) 
plt.colorbar(cax=cbaxes)
plt.savefig("rnaseq_heatmap")
