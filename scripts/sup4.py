import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot
import os

os.system("python ../multimds.py -P 0.1 -w 0 ctrl_Scer_12_32kb.bed galactose_Scer_12_32kb.bed")
struct1 = dt.structure_from_file("ctrl_Scer_12_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Scer_12_32kb_structure.tsv")
plot.plot_structures_interactive((struct1, struct2), out_path="sup4a.png")

os.system("python ../multimds.py -P 0.1 -w 0 ctrl_Scer_12-upstream_32kb.bed galactose_Scer_12-upstream_32kb.bed")
struct1 = dt.structure_from_file("ctrl_Scer_12-upstream_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Scer_12-upstream_32kb_structure.tsv")

plot.plot_structures_interactive((struct1, struct2), out_path="sup4b_upstream.png")

os.system("python ../multimds.py -P 0.1 -w 0 ctrl_Scer_12-downstream_32kb.bed galactose_Scer_12-downstream_32kb.bed")
struct1 = dt.structure_from_file("ctrl_Scer_12-downstream_32kb_structure.tsv")
struct2 = dt.structure_from_file("galactose_Scer_12-downstream_32kb_structure.tsv")

plot.plot_structures_interactive((struct1, struct2), out_path="sup4b_downstream.png")
