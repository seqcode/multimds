set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

PARTITION_NUMS=(4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2)
MIDPOINTS=(135 93 92 51 48 60 60 45 41 53 36 0 0 0 40 24 17 26 28 0 0)
CHROMS=(1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22)
SMOOTHING_PARAMETERS=(6 6 6.5 6 5 5 3 3 3 3 4.5 3 2.5 3 5 4 2.5 3 3 2.5 1)

for i in `seq 0 20`
do
	python ../relocalization_peaks.py GM12878_combined_${CHROMS[$i]}_100kb.bed K562_${CHROMS[$i]}_100kb.bed -m $((${MIDPOINTS[$i]} * 1000000)) -N ${PARTITION_NUMS[$i]} -s ${SMOOTHING_PARAMETERS[$i]} -x ../
done
