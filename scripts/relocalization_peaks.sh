set -e

RES=$1

if [ -e peaks_filtered.bed ]
	then
		rm peaks_filtered.bed
fi

for CHROM in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	python relocalization_peaks.py GM12878_combined K562 $CHROM $RES
	if [ -s ${CHROM}_A_relocalization.bed ]
		then
			./filter_mappability.sh ${CHROM}_A_relocalization
			cat ${CHROM}_A_relocalization_filtered.bed >> peaks_filtered.bed
	fi
done
