set -e

RES=$1
RES_KB=$(($RES/1000))

mkdir -p binding_data

cd binding_data

if [ ! -e GM12878_IDEAS.bb ]
	then
		curl http://bx.psu.edu/~yuzhang/Roadmap_ideas/test1.114.bb -o GM12878_IDEAS.bb
fi

if [ ! -e K562_IDEAS.bb ]
	then
		curl http://bx.psu.edu/~yuzhang/Roadmap_ideas/test1.121.bb -o K562_IDEAS.bb
fi

for CELL_TYPE in GM12878 K562
do
	if [ ! -e ${CELL_TYPE}_IDEAS.bed ]
		then
			bigBedToBed ${CELL_TYPE}_IDEAS.bb ${CELL_TYPE}_IDEAS.bed
	fi	

	if [ ! -e ${CELL_TYPE}_active.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "11_EnhBiv" || $4 == "4_Enh" || $4 == "18_Enh/Het" || $4 == "19_Enh/ReprPC" || $4 == "6_EnhG" || $4 == "17_EnhGA" || $4 == "2_TxWk" || $4 == "5_Tx" || $4 == "14_Tss/Wk" || $4 == "8_TssAFlnk" || $4 == "10_TssA" || $4 == "15_TssBiv" {print $0}' > ${CELL_TYPE}_active.bed
	fi

	WINDOW_FILE=hg19_${RES_KB}kb_windows.bed

	if [ ! -e $WINDOW_FILE ]
		then
			curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -o hg19.chrom.sizes
			bedtools makewindows -g hg19.chrom.sizes -w $RES > $WINDOW_FILE
	fi

	COVERAGE_FILE=$CELL_TYPE"_"${RES_KB}kb_active_coverage.bed

	if [ ! -e $COVERAGE_FILE ]
		then
			bedtools coverage -a $WINDOW_FILE -b $CELL_TYPE"_active".bed > $COVERAGE_FILE
	fi

	for CHROM in `seq 22`
	do
		if [ ! -e $CELL_TYPE"_"$CHROM"_"${RES_KB}kb_active_coverage.bed ]
			then
				cat $COVERAGE_FILE | awk -v chrom=chr$CHROM '$1 == chrom {print $0}' > $CELL_TYPE"_"$CHROM"_"${RES_KB}kb_active_coverage.bed
		fi
	done

done

cd ..
