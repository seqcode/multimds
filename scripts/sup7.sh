set -e

./get_hic_data.sh GM12878_combined 100000
python sup7.py
