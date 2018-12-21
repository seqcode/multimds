set -e

CELL_TYPE=$1
RES=$2
RES_KB=$(($RES/1000))

mkdir -p hic_data

cd hic_data

if [ ! -d $CELL_TYPE/${RES_KB}kb_resolution_intrachromosomal ]
	then
		if [ ! -e GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz ]
			then
				curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz -o GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz
		fi
	tar xzf GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz $CELL_TYPE/${RES_KB}kb_resolution_intrachromosomal
fi

if [ $CELL_TYPE == "K562" ]
	then
		CHROMS=`seq 21`		#skip philadelphia chromosome
	else
		CHROMS=`seq 22`
fi

for CHROM in $CHROMS
do
	echo $CHROM
	if [ -d $CELL_TYPE/${RES_KB}kb_resolution_intrachromosomal/chr$CHROM ] && [ ! -e ${CELL_TYPE}_${CHROM}_${RES_KB}kb.bed ]
		then
			python ../normalize.py $CELL_TYPE $RES $CHROM
	fi
done
	
cd ..
