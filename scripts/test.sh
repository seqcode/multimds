set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

python ../multimds.py --full hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -N 2 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -l 5 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -n 1 hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python ../multimds.py -a 3 hic_data/GM12878_combined_21_10kb.bed hic_data/K562_21_10kb.bed
python ../multimds.py -m 28000000 hic_data/GM12878_combined_20_100kb.bed hic_data/K562_20_100kb.bed
