set -e

RES=$1
RES_KB=$(($RES/1000))

#./get_hic_data.sh GM12878_primary
#./get_hic_data.sh GM12878_replicate
#./get_hic_data.sh K562

#if [ ! -e peaks_filtered.bed ]
#	then 
#		./relocalization_peaks.sh
#fi

python edger_input.py $RES_KB
for CHROM in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	echo $CHROM
	Rscript run_edger.R chr$CHROM_${RES_KB}kb_edgeR_table.tsv chr$CHROM_${RES_KB}kb_edgeR_output.tsv
	python get_sig.py chr$CHROM_${RES_KB}kb_edgeR_output.tsv
done

#unique enhancers
mkdir -p binding_data
cd binding_data

if [ ! -e GM12878_enhancers.bed ] 
	then
		if [ ! -e GM12878.csv ]
			then
				wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012270-mmc7.zip
				unzip 1-s2.0-S0092867413012270-mmc7.zip
		fi
		cat GM12878.csv | awk -F "," '$1 != "track name=\"Enhancers in GM12878\" itemRGB=On color=0" {print $2"\t"$3"\t"$4}' > GM12878_enhancers.bed
fi

if [ ! -e K562_enhancers.bed ] 
	then
		if [ ! -e K562.csv ]
			then
				wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012270-mmc7.zip
				unzip 1-s2.0-S0092867413012270-mmc7.zip
		fi
		cat K562.csv | awk -F "," '$1 != "track name=\"Enhancers in K562\" itemRGB=On color=0" {print $2"\t"$3"\t"$4}' > K562_enhancers.bed
fi

cd ..

if [ ! -e peaks_filtered_GM12878_enhancer_coverage.bed ]
	then
		bedtools coverage -a peaks_filtered.bed -b binding_data/GM12878_enhancers.bed > peaks_filtered_GM12878_enhancer_coverage.bed 
fi

if [ ! -e peaks_filtered_K562_enhancer_coverage.bed ]
	then
		bedtools coverage -a peaks_filtered.bed -b binding_data/K562_enhancers.bed > peaks_filtered_K562_enhancer_coverage.bed 
fi

paste peaks_filtered_GM12878_enhancer_coverage.bed peaks_filtered_K562_enhancer_coverage.bed | awk '$7 > 0.1 && $14 <= 0.1 {print $1"\t"$2"\t"$3}' > peaks_filtered_GM12878_only_enhancer.bed
paste peaks_filtered_GM12878_enhancer_coverage.bed peaks_filtered_K562_enhancer_coverage.bed | awk '$7 <= 0.1 && $14 > 0.1 {print $1"\t"$2"\t"$3}' > peaks_filtered_K562_only_enhancer.bed
paste peaks_filtered_GM12878_enhancer_coverage.bed peaks_filtered_K562_enhancer_coverage.bed | awk '$7 > 0.1 && $14 > 0.1 {print $1"\t"$2"\t"$3}' > peaks_filtered_both_enhancer.bed

#polycomb

cd binding_data

WINDOW_FILE=hg19_${RES_KB}kb_windows.bed
if [ ! -e $WINDOW_FILE ]
	then
		wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
		bedtools makewindows -g hg19.chrom.sizes -w $RES > $WINDOW_FILE
fi

if [ ! -e wgEncodeBroadHistoneGm12878H3k27me3StdPkV2_${RES_KB}kb_windows_enrichment.bed ]
	then 
		wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H3k27me3StdPkV2.broadPeak.gz
		gunzip wgEncodeBroadHistoneGm12878H3k27me3StdPkV2.broadPeak.gz
		
		bedtools coverage -a $WINDOW_FILE -b wgEncodeBroadHistoneGm12878H3k27me3StdPkV2.broadPeak > wgEncodeBroadHistoneGm12878H3k27me3StdPkV2_${RES_KB}kb_windows_enrichment.bed
fi

if [ ! -e wgEncodeBroadHistoneK562H3k27me3StdPk_${RES_KB}kb_windows_enrichment.bed ]
	then 
		wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak.gz
		gunzip wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak.gz
		
		bedtools coverage -a $WINDOW_FILE -b wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak > wgEncodeBroadHistoneK562H3k27me3StdPk_${RES_KB}kb_windows_enrichment.bed
fi

bedtools coverage -a $WINDOW_FILE -b GM12878_enhancers.bed > GM12878_enhancers_${RES_KB}kb_windows_enrichment.bed
bedtools coverage -a $WINDOW_FILE -b K562_enhancers.bed > K562_enhancers_${RES_KB}kb_windows_enrichment.bed

cd ..

#negative control
if [ ! -e A_compartment.bed ]
	then
		python get_a_compartment.py
fi

if [ ! -e A_background.bed ]
	then
		bedtools subtract -a A_compartment.bed -b peaks_filtered.bed > A_background.bed
fi

if [ ! -e A_background_filtered.bed ]
	then
		./filter_mappability.sh A_background
fi

python loop_partners_polycomb.py
