set -e

RES=10000
RES_KB=$(($RES/1000))

#./get_hic_data.sh GM12878_combined $RES
#./get_hic_data.sh K562 $RES
./get_activity_data.sh 100000	#lower resolution (will be used for compartment calculation)
#./relocalization_peaks_independent.sh $RES
./get_state_data.sh

#negative control
if [ ! -e A_compartment_${RES_KB}kb.bed ]
	then
		python get_a_compartment.py $RES
fi

if [ ! -e A_background_independent.bed ]
	then
		bedtools subtract -a A_compartment_${RES_KB}kb.bed -b peaks_filtered_independent.bed | cut -f 1,2,3 > A_background_independent.bed
fi

if [ ! -e A_background_independent_filtered.bed ]
	then
		./filter_mappability.sh A_background_independent $RES
fi

IDS=(H3k27acStdPk H3k04me1StdPkV2 H3k04me3StdPkV2 H3k9acStdPk H3k36me3StdPk H2azStdPk H3k4me2StdPk H3k79me2StdPk H4k20me1StdPk H3k27me3StdPkV2 Ezh239875Pk H3k9me3StdPk H3k36me3StdPk CtcfStdPk)
NAMES=(H3K27ac H3K4me1 H3K4me3 H3K9ac H3K36me3 H2AZ H3K4me2 H3K79me2 H4K20me1 H3K27me3 EZH2 H3K9me3 H3K36me3 CTCF)

for i in `seq 0 $((${#IDS[@]}-1))`
do
	ID=${IDS[$i]}
	NAME=${NAMES[$i]}
	FILENAME=wgEncodeBroadHistoneGm12878${ID}.broadPeak
	if [ ! -e binding_data/$FILENAME ]
		then
			curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/$FILENAME.gz -o binding_data/$FILENAME.gz
			gunzip binding_data/$FILENAME.gz
	fi
	bedtools coverage -a peaks_filtered_independent.bed -b binding_data/$FILENAME > peaks_filtered_independent_GM12878_${NAME}_coverage.bed
	bedtools coverage -a A_background_independent_filtered.bed -b binding_data/$FILENAME > A_background_independent_filtered_GM12878_${NAME}_coverage.bed
done

IDS=(H3k27acStdPk H3k4me1StdPk H3k4me3StdPk H3k9acStdPk H3k36me3StdPk H2azStdPk H3k4me2StdPk H3k79me2StdPk H4k20me1StdPk H3k27me3StdPk Ezh239875StdPk H3k9me3StdPk H3k36me3StdPk CtcfStdPk)
NAMES=(H3K27ac H3K4me1 H3K4me3 H3K9ac H3K36me3 H2AZ H3K4me2 H3K79me2 H4K20me1 H3K27me3 EZH2 H3K9me3 H3K36me3 CTCF)

for i in `seq 0 $((${#IDS[@]}-1))`
do
	ID=${IDS[$i]}
	NAME=${NAMES[$i]}
	FILENAME=wgEncodeBroadHistoneK562${ID}.broadPeak
	if [ ! -e binding_data/$FILENAME ]
		then
			curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/$FILENAME.gz -o binding_data/$FILENAME.gz
			gunzip binding_data/$FILENAME.gz
	fi
	bedtools coverage -a peaks_filtered_independent.bed -b binding_data/$FILENAME > peaks_filtered_independent_K562_${NAME}_coverage.bed
	bedtools coverage -a A_background_independent_filtered.bed -b binding_data/$FILENAME > A_background_independent_filtered_K562_${NAME}_coverage.bed
done

for STATE in promoter poised_promoter enhancer insulator transcription repressed heterochromatin
do
	for CELLTYPE in Gm12878 K562
	do
		bedtools coverage -a peaks_filtered_independent.bed -b binding_data/${CELLTYPE}_${STATE}.bed > peaks_filtered_independent_${CELLTYPE}_${STATE}_coverage.bed
		bedtools coverage -a A_background_independent_filtered.bed -b binding_data/${CELLTYPE}_${STATE}.bed > A_background_independent_filtered_${CELLTYPE}_${STATE}_coverage.bed
	done
done

python relocalization_enrichment_independent.py
