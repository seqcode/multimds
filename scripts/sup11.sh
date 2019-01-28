set -e

RES=100000
./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES

python compartment_relocalization_line_plot.py
