set -e

RES=$1
RES_KB=$(($RES/1000))

PARTITION_NUMS=(4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2)
MIDPOINTS=(135 93 92 51 48 60 60 45 41 53 36 0 0 0 40 24 17 26 28 0 0)
CHROMS=(1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 X)
SMOOTHING_PARAMETERS=(6 6 6.5 6 5 5 3 3 3 3 4.5 3 2.5 3 5 4 2.5 3 3 2.5 1)
#SMOOTHING_PARAMETERS=(3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3)

if [ -e peaks_filtered.bed ]
	then
		rm peaks_filtered.bed
fi

for i in `seq 0 19`
do
	CHROM=${CHROMS[$i]}
	#python relocalization_peaks.py GM12878_combined K562 $CHROM $((${MIDPOINTS[$i]} * 1000000)) ${PARTITION_NUMS[$i]} ${SMOOTHING_PARAMETERS[$i]} $RES
	#bedtools subtract -A -a ${CHROM}_dist_peaks.bed -b ${CHROM}_comp_peaks.bed > ${CHROM}_noncomp_peaks.bed
	cat ${CHROM}_dist_peaks.bed | awk '($4 - $5 < 0.2 && $4 > $5) || ($5 - $4 < 0.2 && $5 > $4) {print $0}' > ${CHROM}_noncomp_peaks.bed
	#cat ${CHROM}_noncomp_peaks.bed | awk '$4 > 0 && $5 > 0 {print $1"\t"$2"\t"$3}' > ${CHROM}_A_noncomp_peaks.bed	#A compartment only
	cat ${CHROM}_noncomp_peaks.bed | awk '$4 > 0 && $5 > 0 {print $0}' > ${CHROM}_A_noncomp_peaks.bed	#A compartment only
	./filter_mappability.sh ${CHROM}_A_noncomp_peaks $RES
	cat ${CHROM}_A_noncomp_peaks_filtered.bed >> peaks_filtered.bed
done
