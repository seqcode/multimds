set -e

./get_hic_data.sh GM12878_combined 10000
./get_hic_data.sh K562 10000

python sup24.py
