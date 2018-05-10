set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

python ../multimds.py hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed -N 2
