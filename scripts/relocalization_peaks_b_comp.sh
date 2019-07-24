set -e

RES=$1
RES_KB=$(($RES/1000))

#if [ -e peaks_b_comp_filtered.bed ]
#	then
#		rm peaks_b_comp_filtered.bed
#fi

#for CHROM in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21
#do
#	python ../multimds.py hic_data/GM12878_combined_${CHROM}_${RES_KB}kb.bed hic_data/K562_${CHROM}_${RES_KB}kb.bed
	#python relocalization_peaks.py GM12878_combined K562 $CHROM $RES
	#cat ${CHROM}_dist_peaks.bed | awk '($4 - $5 < 0.2 && $4 > $5) || ($5 - $4 < 0.2 && $5 > $4) {print $0}' > ${CHROM}_noncomp_peaks.bed
#	cat ${CHROM}_noncomp_peaks.bed | awk '$4 < 0 && $5 < 0 {print $0}' > ${CHROM}_B_noncomp_peaks.bed	#A compartment only
#	./filter_mappability.sh ${CHROM}_B_noncomp_peaks $RES
#	cat ${CHROM}_B_noncomp_peaks_filtered.bed >> peaks_b_comp_filtered.bed
#done

#negative control
if [ -e B_compartment_${RES_KB}kb.bed ]
	then
		rm B_compartment_${RES_KB}kb.bed
fi

python get_b_compartment.py $RES
bedtools subtract -a B_compartment_${RES_KB}kb.bed -b peaks_b_comp_filtered.bed > B_background.bed
./filter_mappability.sh B_background $RES
