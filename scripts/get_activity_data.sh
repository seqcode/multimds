set -e

RES=$1
RES_KB=$(($RES/1000))

mkdir -p binding_data

cd binding_data

for CELL_TYPE in Gm12878 K562
do

	if [ ! -e $CELL_TYPE"_"*_${RES}kb_active_coverage.bed ]
		then
			curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz -o wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz
			gunzip wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "1_Active_Promoter" || $4 == "2_Weak_Promoter" || $4 == "3_Poised_Promoter" || $4 == "4_Strong_Enhancer" || $4 == "5_Strong_Enhancer" || $4 == "6_Weak_Enhancer" || $4 == "7_Weak_Enhancer" || $4 == "9_Txn_Transition" || $4 == "10_Txn_Elongation" || $4 == "11_Weak_Txn" {print $0}' > $CELL_TYPE"_active".bed

		WINDOW_FILE=hg19_${RES_KB}kb_windows.bed

		if [ ! -e $WINDOW_FILE ]
			then
				curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -o hg19.chrom.sizes
				bedtools makewindows -g hg19.chrom.sizes -w $RES > $WINDOW_FILE
		fi

		COVERAGE_FILE=$CELL_TYPE"_"${RES_KB}kb_active_coverage.bed

		bedtools coverage -a $WINDOW_FILE -b $CELL_TYPE"_active".bed > $COVERAGE_FILE

		for CHROM in `seq 22`
		do
			cat $COVERAGE_FILE | awk -v chrom=chr$CHROM '$1 == chrom {print $0}' > $CELL_TYPE"_"$CHROM"_"${RES_KB}kb_active_coverage.bed
		done
	fi
done

cd ..
