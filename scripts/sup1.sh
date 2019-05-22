set -e

RES=100000

./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES

python embedding_error.py

mv GM12878_combined_K562_embedding_error.png sup1.png
