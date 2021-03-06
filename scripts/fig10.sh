set -e

RES=10000
RES_KB=$(($RES/1000))

./get_hic_data.sh GM12878_combined $RES
./get_hic_data.sh K562 $RES
./get_activity_data.sh 100000	#lower resolution (will be used for compartment calculation)
./relocalization_peaks_a_comp.sh $RES

IDS=(H3k27acStdPk H3k04me1StdPkV2 H3k04me3StdPkV2 H3k9acStdPk H3k36me3StdPk H2azStdPk H3k4me2StdPk H3k79me2StdPk H4k20me1StdPk H3k27me3StdPkV2 Ezh239875Pk H3k9me3StdPk H3k36me3StdPk CtcfStdPk)
NAMES=(H3K27ac H3K4me1 H3K4me3 H3K9ac H3K36me3 H2AZ H3K4me2 H3K79me2 H4K20me1 H3K27me3 EZH2 H3K9me3 H3K36me3 CTCF)

for i in `seq 0 $((${#IDS[@]}-1))`
do
	ID=${IDS[$i]}
	NAME=${NAMES[$i]}
	FILENAME=wgEncodeBroadHistoneGm12878${ID}.broadPeak
	if [ ! -s binding_data/$FILENAME ]
		then
			curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/$FILENAME.gz -o binding_data/$FILENAME.gz
			gunzip binding_data/$FILENAME.gz
	fi
	bedtools coverage -a peaks_A_filtered.bed -b binding_data/$FILENAME > peaks_A_filtered_GM12878_${NAME}_coverage.bed
done

IDS=(H3k27acStdPk H3k4me1StdPk H3k4me3StdPk H3k9acStdPk H3k36me3StdPk H2azStdPk H3k4me2StdPk H3k79me2StdPk H4k20me1StdPk H3k27me3StdPk Ezh239875StdPk H3k9me3StdPk H3k36me3StdPk CtcfStdPk)
NAMES=(H3K27ac H3K4me1 H3K4me3 H3K9ac H3K36me3 H2AZ H3K4me2 H3K79me2 H4K20me1 H3K27me3 EZH2 H3K9me3 H3K36me3 CTCF)

for i in `seq 0 $((${#IDS[@]}-1))`
do
	ID=${IDS[$i]}
	NAME=${NAMES[$i]}
	FILENAME=wgEncodeBroadHistoneK562${ID}.broadPeak
	if [ ! -s binding_data/$FILENAME ]
		then
			curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/$FILENAME.gz -o binding_data/$FILENAME.gz
			gunzip binding_data/$FILENAME.gz
	fi
	bedtools coverage -a peaks_A_filtered.bed -b binding_data/$FILENAME > peaks_A_filtered_K562_${NAME}_coverage.bed
done

paste peaks_A_filtered_GM12878_H3K27ac_coverage.bed peaks_A_filtered_GM12878_H3K4me1_coverage.bed peaks_A_filtered_GM12878_H3K4me3_coverage.bed peaks_A_filtered_GM12878_H3K9ac_coverage.bed peaks_A_filtered_GM12878_H2AZ_coverage.bed peaks_A_filtered_GM12878_H2AZ_coverage.bed peaks_A_filtered_GM12878_H3K4me2_coverage.bed peaks_A_filtered_GM12878_H3K27me3_coverage.bed peaks_A_filtered_K562_H3K27ac_coverage.bed peaks_A_filtered_K562_H3K4me1_coverage.bed peaks_A_filtered_K562_H3K4me3_coverage.bed peaks_A_filtered_K562_H3K9ac_coverage.bed peaks_A_filtered_K562_H2AZ_coverage.bed peaks_A_filtered_K562_H3K4me2_coverage.bed peaks_A_filtered_K562_H3K27me3_coverage.bed | cut -f 9,18,27,36,45,54,63,72,81,90,99,108,117,126,135 > peaks_filtered_all_coverage.bed

python fig10.py
