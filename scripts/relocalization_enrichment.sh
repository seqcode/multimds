set -e

RES=10000
RES_KB=$(($RES/1000))

#./get_hic_data.sh GM12878_combined $RES
#./get_hic_data.sh K562 $RES
#./get_activity_data.sh
#./relocalization_peaks.sh $RES

#negative control
if [ ! -e A_compartment_${RES_KB}kb.bed ]
	then
		python get_a_compartment.py $RES
fi

if [ ! -e A_background.bed ]
	then
		bedtools subtract -a A_compartment_${RES_KB}kb.bed -b peaks_filtered.bed > A_background.bed
fi

if [ ! -e A_background_filtered.bed ]
	then
		./filter_mappability.sh A_background $RES
fi

IDS=(H3k27acStdPk H3k04me1StdPkV2 H3k04me3StdPkV2 H3k9acStdPk CtcfStdPk H2azStdPk H3k4me2StdPk H3k79me2StdPk H4k20me1StdPk)
NAMES=(H3K27ac H3K4me1 H3K4me3 H3K9ac CTCF H2AZ H3K4me2 H3K79me2 H4K20me1)

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
	bedtools coverage -a peaks_filtered.bed -b binding_data/$FILENAME > peaks_filtered_GM12878_${NAME}_coverage.bed
	bedtools coverage -a A_background_filtered.bed -b binding_data/$FILENAME > A_background_filtered_GM12878_${NAME}_coverage.bed
	#python chi_square.py $(cat peaks_filtered_GM12878_${NAME}_coverage.bed | awk '$7 > 0 {print 1}' | wc -l) $(cat peaks_filtered_GM12878_${NAME}_coverage.bed | awk '$7 == 0 {print 1}' | wc -l) $(cat A_background_filtered_GM12878_${NAME}_coverage.bed | awk '$7 > 0 {print 1}' | wc -l) $(cat A_background_filtered_GM12878_${NAME}_coverage.bed | awk '$7 == 0 {print 1}' | wc -l) GM12878 $NAME
done

NAME=H3K27me3
FILENAME=wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak
if [ ! -e binding_data/$FILENAME ]
	then
		curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/$FILENAME.gz -o binding_data/$FILENAME.gz
		gunzip binding_data/$FILENAME.gz
fi
bedtools coverage -a peaks_filtered.bed -b binding_data/$FILENAME > peaks_filtered_K562_${NAME}_coverage.bed
bedtools coverage -a A_background_filtered.bed -b binding_data/$FILENAME > A_background_filtered_K562_${NAME}_coverage.bed
#python chi_square.py $(cat peaks_filtered_K562_${NAME}_coverage.bed | awk '$7 > 0 {print 1}' | wc -l) $(cat peaks_filtered_K562_${NAME}_coverage.bed | awk '$7 == 0 {print 1}' | wc -l) $(cat A_background_filtered_K562_${NAME}_coverage.bed | awk '$7 > 0 {print 1}' | wc -l) $(cat A_background_filtered_K562_${NAME}_coverage.bed | awk '$7 == 0 {print 1}' | wc -l) K562 $NAME

python chi_square.py
