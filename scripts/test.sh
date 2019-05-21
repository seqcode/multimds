set -e

PYTHON=$1

./get_hic_data.sh GM12878_combined 100000
./get_hic_data.sh K562 100000

$PYTHON ../multimds.py hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py -P 0.5 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py -w 0.5 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py --partitioned hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py --partitioned -N 4 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py --partitioned -l 5 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py -n 1 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON ../multimds.py -a 3 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON test_plot.py
rm *structure.tsv
$PYTHON ../multimds.py -o test_ hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
$PYTHON ../multimds.py --partitioned -r 28000000 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
