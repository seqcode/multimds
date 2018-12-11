set -e

RES=$1

#./get_hic_data.sh GM12878_combined $RES
#./get_hic_data.sh K562 $RES
#./get_activity_data.sh $RES
./relocalization_peaks.sh $RES

#enhancers
if [ ! -e GM12878_enhancers.bed ] 
	then
		if [ ! -e GM12878.csv ]
			then
				wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012270-mmc7.zip
				unzip 1-s2.0-S0092867413012270-mmc7.zip
		fi
		cat GM12878.csv | awk -F "," '$1 != "track name=\"Enhancers in GM12878\" itemRGB=On color=0" {print $2"\t"$3"\t"$4}' > GM12878_enhancers.bed
fi

#if [ ! -e peaks_filtered_GM12878_enhancer_coverage.bed ]
#	then
		bedtools coverage -a peaks_filtered.bed -b GM12878_enhancers.bed > peaks_filtered_GM12878_enhancer_coverage.bed 
#fi

#negative control
if [ ! -e A_compartment.bed ]
	then
		python get_a_compartment.py
fi

#if [ ! -e A_background.bed ]
#	then
		bedtools subtract -a A_compartment.bed -b peaks_filtered.bed > A_background.bed
#fi

#if [ ! -e A_background_filtered.bed ]
#	then
		./filter_mappability.sh A_background
#fi

#if [ ! -e A_background_filtered_GM12878_enhancer_coverage.bed ]
#	then
		bedtools coverage -a A_background_filtered.bed -b GM12878_enhancers.bed > A_background_filtered_GM12878_enhancer_coverage.bed 
#fi

python enhancer_pie.py $(cat peaks_filtered_GM12878_enhancer_coverage.bed | awk '$7 > 0.05 {print 1}' | wc -l) $(cat peaks_filtered_GM12878_enhancer_coverage.bed | awk '$7 <= 0.05 {print 1}' | wc -l) $(cat A_background_filtered_GM12878_enhancer_coverage.bed | awk '$7 > 0.05 {print 1}' | wc -l) $(cat A_background_filtered_GM12878_enhancer_coverage.bed | awk '$7 <= 0.05 {print 1}' | wc -l)

python ttest.py peaks_filtered_GM12878_enhancer_coverage.bed A_background_filtered_GM12878_enhancer_coverage.bed
