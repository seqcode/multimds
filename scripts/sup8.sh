set -e

./get_hic_data.sh GM12878_combined 100000
python ../minimds.py hic_data/GM12878_combined_21_100kb.bed
python sup8.py
