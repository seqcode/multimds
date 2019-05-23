set -e

RES=10000
./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES
./relocalization_peaks.sh $RES

python peaks_vs_background_compartments.py
