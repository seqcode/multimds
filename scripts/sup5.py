import numpy as np
from multimds import plotting as plot
from multimds import multimds as mm

mappability = np.loadtxt("mappability_21_5kb.bed", usecols=3)

#struct1 = dt.structure_from_file("GM12878_combined_21_5kb_structure.tsv")
#struct2 = dt.structure_from_file("K562_21_5kb_structure.tsv")

struct1, struct2 = mm.full_mds("hic_data/GM12878_combined_21_5kb.bed", "hic_data/K562_21_5kb.bed")

mappability = mappability[struct1.chrom.minPos/struct1.chrom.res + struct1.nonzero_abs_indices()]	#only loci in structures
mappable = np.where(mappability > 0.75)
struct1.points = struct1.points[mappable]
struct2.points = struct2.points[mappable]

plot.plot_structures_interactive((struct1, struct2))
