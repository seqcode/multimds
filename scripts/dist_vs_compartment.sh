set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562
#./get_activity_data.sh

python dist_vs_compartment.py
