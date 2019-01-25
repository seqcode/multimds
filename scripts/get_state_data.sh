set -e

mkdir -p binding_data

cd binding_data

for CELL_TYPE in Gm12878 K562
do
	if [ ! -e wgEncodeBroadHmm$CELL_TYPE"HMM".bed ]
		then
			curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz -o wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz
			gunzip wgEncodeBroadHmm$CELL_TYPE"HMM".bed.gz
	fi

	if [ ! -e $CELL_TYPE"_promoter".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "1_Active_Promoter" || $4 == "2_Weak_Promoter" {print $0}' > $CELL_TYPE"_promoter".bed
	fi

	if [ ! -e $CELL_TYPE"_poised_promoter".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "3_Poised_Promoter" {print $0}' > $CELL_TYPE"_poised_promoter".bed
	fi

	if [ ! -e $CELL_TYPE"_enhancer".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "4_Strong_Enhancer" || $4 == "5_Strong_Enhancer" || $4 == "6_Weak_Enhancer" || $4 == "7_Weak_Enhancer" {print $0}' > $CELL_TYPE"_enhancer".bed
	fi

	if [ ! -e $CELL_TYPE"_insulator".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "8_Insulator" {print $0}' > $CELL_TYPE"_insulator".bed
	fi

	if [ ! -e $CELL_TYPE"_transcription".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "9_Txn_Transition" || $4 == "10_Txn_Elongation" || $4 == "11_Weak_Txn" {print $0}' > $CELL_TYPE"_transcription".bed
	fi

	if [ ! -e $CELL_TYPE"_repressed".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "12_Repressed" {print $0}' > $CELL_TYPE"_repressed".bed
	fi

	if [ ! -e $CELL_TYPE"_heterochromatin".bed ]
		then
			cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "13_Heterochrom/lo" || $4 == "14_Repetitive/CNV" || $4 == "15_Repetitive/CNV" {print $0}' > $CELL_TYPE"_heterochromatin".bed
	fi
	

done

cd ..
