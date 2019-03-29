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
done


if [ ! -e Gm12878_enhancer.bed ]
	then
		cat wgEncodeBroadHmmGm12878HMM.bed | awk '$4 == "4_Strong_Enhancer" || $4 == "5_Strong_Enhancer" || $4 == "6_Weak_Enhancer" || $4 == "7_Weak_Enhancer" {print $0}' > Gm12878_enhancer.bed
fi


if [ ! -e $CELL_TYPE"_repressed".bed ]
	then
		cat wgEncodeBroadHmm$CELL_TYPE"HMM".bed | awk '$4 == "12_Repressed" {print $0}' > $CELL_TYPE"_repressed".bed
fi

cd ..
