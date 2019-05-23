set -e

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

	if [ ! -e ${CELL_TYPE}_enhancer.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "11_EnhBiv" || $4 == "4_Enh" || $4 == "18_Enh/Het" || $4 == "19_Enh/ReprPC" || $4 == "6_EnhG" || $4 == "17_EnhGA" {print $0}' > ${CELL_TYPE}_enhancer.bed
	fi


	if [ ! -e ${CELL_TYPE}_polycomb.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "1_ReprPCWk" || $4 == "12_Het/ReprPC" || $4 == "13_ReprPC" {print $0}' > ${CELL_TYPE}_polycomb.bed
	fi

	if [ ! -e ${CELL_TYPE}_transcription.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "2_TxWk" || $4 == "5_Tx" {print $0}' > ${CELL_TYPE}_transcription.bed
	fi
	if [ ! -e ${CELL_TYPE}_TSS.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "14_Tss/Wk" || $4 == "8_TssAFlnk" || $4 == "10_TssA" || $4 == "15_TssBiv" {print $0}' > ${CELL_TYPE}_TSS.bed
	fi


done

cd ..
