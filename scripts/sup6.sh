set -e

./get_yeast_data.sh

python ../multimds.py -P 0.1 -w 0 ctrl_Scer_12_32kb.bed galactose_Scer_12_32kb.bed
python ../multimds.py -P 0.1 -w 0 ctrl_Scer_12-upstream_32kb.bed galactose_Scer_12-upstream_32kb.bed
python ../multimds.py -P 0.1 -w 0 ctrl_Scer_12-downstream_32kb.bed galactose_Scer_12-downstream_32kb.bed
python sup6.py
