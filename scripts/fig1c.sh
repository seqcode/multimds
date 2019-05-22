set -e

RES=100000

./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES

python reproducibility.py GM12878_combined_21_100kb.bed K562_21_100kb.bed
mv GM12878_combined_K562_reproducibility.png fig1c.png
