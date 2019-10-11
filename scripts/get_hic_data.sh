set -e

CELL_TYPE=$1
RES=$2
RES_KB=$(($RES/1000))

mkdir -p hic_data

cd hic_data

if [ ! -d $CELL_TYPE/${RES_KB}kb_resolution_intrachromosomal ]
	then
		if [ ! -s GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz ]
			then
				curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz -o GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz
		fi
	tar xzf GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz $CELL_TYPE/${RES_KB}kb_resolution_intrachromosomal
fi

if [ $CELL_TYPE == "K562" ] || [ $CELL_TYPE == "KBM7" ]
	then
		CHROMS=(1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21)	#skip translocated
	else
		CHROMS=`seq 22`
fi

for CHROM in ${CHROMS[*]}
do
	if [ -d $CELL_TYPE/${RES_KB}kb_resolution_intrachromosomal/chr$CHROM ] && [ ! -s ${CELL_TYPE}_${CHROM}_${RES_KB}kb.bed ]
		then
			echo $CHROM
			python ../normalize.py $CELL_TYPE $RES $CHROM
	fi
done
	
cd ..
