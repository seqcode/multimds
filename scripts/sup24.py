from multimds import plotting as plot
from multimds import multimds as mm

struct1, struct2 = mm.full_mds("hic_data/GM12878_combined_21_10kb.bed", "hic_data/K562_21_10kb.bed", weight=0)

plot.plot_structures_interactive((struct1, struct2))
