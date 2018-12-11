from matplotlib import pyplot as plt
import sys
sys.path.append("..")
import compartment_analysis as ca
import data_tools as dt
import os

paths = sys.argv[1:len(sys.argv)]
prefixes = [os.path.basename(path) for path in paths]
structs = [dt.structureFromBed(path) for path in paths]
mats = [dt.matFromBed(path, struct) for path, struct in zip(paths, structs)]
all_comps = [ca.get_compartments(mat) for mat in mats]
all_gen_coords = [struct.getGenCoords() for struct in structs]

#all_comps[len(all_comps)-1] = -all_comps[len(all_comps)-1]

for gen_coords, comps, prefix in zip(all_gen_coords, all_comps, prefixes):
	plt.plot(gen_coords, comps, label=prefix)

plt.legend()
plt.show()
