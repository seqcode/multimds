set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

python embedding_error.py 21
