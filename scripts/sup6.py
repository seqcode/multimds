import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot

struct1 = dt.structure_from_file("ctrl_Scer_12_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Scer_12_32kb_structure.tsv")
plot.plot_structures_interactive((struct1, struct2), out_path="sup6a.png")

struct1 = dt.structure_from_file("ctrl_Scer_12-upstream_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Scer_12-upstream_32kb_structure.tsv")
plot.plot_structures_interactive((struct1, struct2), out_path="sup6b_upstream.png")

struct1 = dt.structure_from_file("ctrl_Scer_12-downstream_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Scer_12-downstream_32kb_structure.tsv")
plot.plot_structures_interactive((struct1, struct2), out_path="sup6b_downstream.png")
