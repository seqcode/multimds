from multimds import plotting as plot
from multimds import multimds as mm

struct1, struct2 = mm.full_mds("ctrl_Scer_12_32kb.bed", "galactose_Scer_12_32kb.bed", weight=0, penalty=0.1)
plot.plot_structures_interactive((struct1, struct2), out_path="sup6a.png")

struct1, struct2 = mm.full_mds("ctrl_Scer_12-upstream_32kb.bed", "galactose_Scer_12-upstream_32kb.bed", weight=0, penalty=0.1)
plot.plot_structures_interactive((struct1, struct2), out_path="sup6b_upstream.png")

struct1, struct2 = mm.full_mds("ctrl_Scer_12-downstream_32kb.bed", "galactose_Scer_12-downstream_32kb.bed", weight=0, penalty=0.1)
plot.plot_structures_interactive((struct1, struct2), out_path="sup6b_downstream.png")
