set -e

RES=$1
RES_KB=$(($RES/1000))

#./get_hic_data.sh K562_combined $RES
#./get_hic_data.sh K562 $RES
#./get_activity_data.sh
#./relocalization_peaks.sh $RES

#enhancers
#if [ ! -e K562_enhancers.bed ] 
#	then
#		if [ ! -e K562.csv ]
#			then
#				wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012270-mmc7.zip
#				unzip 1-s2.0-S0092867413012270-mmc7.zip
#		fi
#		cat K562.csv | awk -F "," '$1 != "track name=\"Enhancers in K562\" itemRGB=On color=0" {print $2"\t"$3"\t"$4}' > K562_enhancers.bed
#fi

if [ ! -e peaks_filtered_K562_H3K27me3_coverage.bed ]
	then
		bedtools coverage -a peaks_filtered.bed -b binding_data/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak > peaks_filtered_K562_H3K27me3_coverage.bed
fi

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

if [ ! -e A_background_filtered_K562_H3K27me3_coverage.bed ]
	then
		bedtools coverage -a A_background_filtered.bed -b binding_data/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak > A_background_filtered_K562_H3K27me3_coverage.bed
fi

python ttest.py peaks_filtered_K562_H3K27me3_coverage.bed A_background_filtered_K562_H3K27me3_coverage.bed K562 H3K27me3
