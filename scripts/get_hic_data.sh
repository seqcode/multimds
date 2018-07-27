set -e

MINIMDS_DIR=$1
CELL_TYPE=$2

mkdir -p hic_data

cd hic_data

if [ ! -e $CELL_TYPE"_"*_100kb.bed ]
	then
		if [ ! -e GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz ]
			then
				wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz
		fi

		tar xzf GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz $CELL_TYPE/100kb_resolution_intrachromosomal
	
		for CHROM in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22	#skip chr9
		do
			echo $CHROM
			python $MINIMDS_DIR/scripts/normalize.py $CELL_TYPE 100000 $CHROM
		done
			
		rm GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz

		cat $CELL_TYPE"_22_100kb".bed | awk '$2 >= 24000000 && $5 >= 24000000 {print $0}' > tmp	#filter out translocation
		mv tmp $CELL_TYPE"_22_100kb".bed
fi	

cd ..
