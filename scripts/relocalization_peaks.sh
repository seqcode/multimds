set -e

RES=100000

#./get_hic_data.sh GM12878_combined
#./get_hic_data.sh K562
#./get_activity_data.sh

PARTITION_NUMS=(4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2)
MIDPOINTS=(135 93 92 51 48 60 60 45 41 53 36 0 0 0 40 24 17 26 28 0 0)
CHROMS=(1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22)
SMOOTHING_PARAMETERS=(6 6 6.5 6 5 5 3 3 3 3 4.5 3 2.5 3 5 4 2.5 3 3 2.5 2)

if [ -e peaks_filtered.bed ]
	then
		rm peaks_filtered.bed
fi

for i in `seq 0 20`
do
	CHROM=${CHROMS[$i]}
	if [ $CHROM -eq 22 ]
		then
			CELL_TYPE1=GM12878_combined_filtered
		else
			CELL_TYPE1=GM12878_combined
	fi
	python relocalization_peaks.py $CELL_TYPE1 K562 $CHROM $((${MIDPOINTS[$i]} * 1000000)) ${PARTITION_NUMS[$i]} ${SMOOTHING_PARAMETERS[$i]}
	bedtools subtract -A -a ${CHROM}_dist_peaks.bed -b ${CHROM}_comp_peaks.bed > ${CHROM}_noncomp_peaks.bed
	cat ${CHROM}_noncomp_peaks.bed | awk -v res=$RES '$5 > 0 && $6 > 0 {print $1"\t"$4"\t"$4+res}' > ${CHROM}_A_noncomp_peaks.bed	#A compartment only
	./filter_mappability.sh ${CHROM}_A_noncomp_peaks
	cat ${CHROM}_A_noncomp_peaks_filtered.bed >> peaks_filtered.bed
done
