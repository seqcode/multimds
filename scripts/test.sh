set -e

./get_hic_data.sh GM12878_combined 100000
./get_hic_data.sh K562 100000

python ../multimds.py hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -P 0.1 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -w 0.1 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py --partitioned hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py --partitioned -N 2 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py --partitioned -l 5 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -n 1 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -a 3 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -o test_ hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py --partitioned -r 28000000 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
