set -e

mkdir -p binding_data

cd binding_data

for CELL_TYPE in Gm12878 K562
do

	if [ ! -s $CELL_TYPE"_"*_100kb_active_coverage.bed ]
		then
			wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz
			gunzip wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "1_Active_Promoter" || $4 == "2_Weak_Promoter" || $4 == "3_Poised_Promoter" || $4 == "4_Strong_Enhancer" || $4 == "5_Strong_Enhancer" || $4 == "6_Weak_Enhancer" || $4 == "7_Weak_Enhancer" || $4 == "9_Txn_Transition" || $4 == "10_Txn_Elongation" || $4 == "11_Weak_Txn" {print $0}' > $CELL_TYPE"_active".bed

		WINDOW_FILE=hg19_100kb_windows.bed

		if [ ! -s $WINDOW_FILE ]
			then
				wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
				bedtools makewindows -g hg19.chrom.sizes -w 100000 > $WINDOW_FILE
		fi

		bedtools coverage -a $WINDOW_FILE -b $CELL_TYPE"_active".bed > $COVERAGE_FILE

		for CHROM in `seq 22`
		do
			cat $COVERAGE_FILE | awk -v chrom=chr$CHROM '$1 == chrom {print $0}' > $CELL_TYPE"_"$CHROM"_100kb_active_coverage".bed
		done
	fi
done

cd ..
