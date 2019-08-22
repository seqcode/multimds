set -e

RES=10000
./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES

./relocalization_peaks_a_comp.sh $RES
python sup21a.py

./relocalization_peaks_b_comp.sh $RES
python sup21b.py
