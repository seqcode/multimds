""""Convert fixedStep wig to binned bed"""

import sys
sys.path.append("..")
from tools import Tracker

wig = sys.argv[1]
bin_size = int(sys.argv[2])
file_size = int(sys.argv[3])

prefix = wig.split(".")[0]

tracker = Tracker("Converting {}".format(wig), file_size)

tot = 0
count = 0

with open(wig) as in_file:
	with open("{}.bed".format(prefix), "w") as out_file:
		for line in in_file:
			line = line.strip().split()
			if line[0] == "fixedStep":	#header
				chrom = line[1].split("=")[1]
				curr_pos = int(line[2].split("=")[1])
				step = int(line[3].split("=")[1])
				span = int(line[4].split("=")[1])
			else:
				tot += float(line[0])
				count += span
				if curr_pos%bin_size == 0:
					if count == 0:
						avg = 0
					else:
						avg = tot/count
					out_file.write("\t".join((chrom, str(curr_pos-bin_size), str(curr_pos), str(avg))))
					out_file.write("\n")
					tot = 0		#re-initialize
					count = 0
				curr_pos += step
			tracker.increment()
		out_file.close()
	in_file.close()
