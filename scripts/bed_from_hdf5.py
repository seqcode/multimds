import h5py
import sys
sys.path.append("/home/lur159/git/miniMDS")
from tools import Tracker
import numpy as np

in_path = sys.argv[1]
out_path = sys.argv[2]
f = h5py.File(in_path)

counts = np.array(f["pixels"]["count"])
bin_ids1 = np.array(f["pixels"]["bin1_id"])
bin_ids2 = np.array(f["pixels"]["bin2_id"])
chroms = np.array(f["bins"]["chrom"])
starts = np.array(f["bins"]["start"])
ends = np.array(f["bins"]["end"])

f.close()

tracker = Tracker("Converting to BED", len(counts))

print "Begin converting to BED"
with open(out_path, "w") as out_file:
	for count, bin_id1, bin_id2 in zip(counts, bin_ids1, bin_ids2):
		if count != 0:
			chrom1 = str(chroms[bin_id1] + 1)	#switch to 1-indexed
			chrom2 = str(chroms[bin_id2] + 1)
			start1 = str(starts[bin_id1])
			end1 = str(ends[bin_id1])
			start2 = str(starts[bin_id2])
			end2 = str(ends[bin_id2])
			out_file.write("\t".join(("chr" + chrom1, start1, end1, "chr" + chrom2, start2, end2, str(count))))
			out_file.write("\n")
		tracker.increment()	
	out_file.close()
