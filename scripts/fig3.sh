set -e

RES=100000
./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES
./get_activity_data.sh $RES

python dist_vs_compartment.py
