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


	if [ ! -e ${CELL_TYPE}_polycomb.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "1_ReprPCWk" || $4 == "12_Het/ReprPC" || $4 == "13_ReprPC" {print $0}' > ${CELL_TYPE}_polycomb.bed
	fi

	if [ ! -e ${CELL_TYPE}_heterochromatin.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "3_HetWk" {print $0}' > ${CELL_TYPE}_heterochromatin.bed
	fi

	if [ ! -e ${CELL_TYPE}_repeats.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "7_ZNF/Rpts" {print $0}' > ${CELL_TYPE}_repeats.bed
	fi

	if [ ! -e ${CELL_TYPE}_quiescent.bed ]
		then
			cat ${CELL_TYPE}_IDEAS.bed | awk '$4 == "0_Quies" {print $0}' > ${CELL_TYPE}_quiescent.bed
	fi

done

cd ..
