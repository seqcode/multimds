set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

python reproducibility.py GM12878_combined_21_100kb.bed K562_21_100kb.bed
