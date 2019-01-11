import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot

struct1 = dt.structure_from_file("GM12878_combined_21_100kb_structure.tsv")
struct2 = dt.structure_from_file("K562_21_100kb_structure.tsv")

plot.plot_structures_interactive((struct1, struct2))
