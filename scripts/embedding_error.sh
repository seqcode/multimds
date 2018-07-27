set -e

MINIMDS_DIR=$1

./get_hic_data.sh $MINIMDS_DIR GM12878_combined
./get_hic_data.sh $MINIMDS_DIR K562

python embedding_error.py 21
