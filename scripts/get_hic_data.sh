set -e

CELL_TYPE=$1

mkdir -p hic_data

cd hic_data

if [ ! -e $CELL_TYPE"_"*_100kb.bed ]
	then
		if [ ! -e GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz ]
			then
				wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz
		fi
		tar xzf GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz $CELL_TYPE/100kb_resolution_intrachromosomal
	
		for CHROM in `seq 22`
		do
			echo $CHROM
			if [ -d $CELL_TYPE/100kb_resolution_intrachromosomal/chr$CHROM ]
				then
					python ../normalize.py $CELL_TYPE 100000 $CHROM
			fi
		done
			
		rm GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz

		cat $CELL_TYPE"_22_100kb".bed | awk '$2 >= 24000000 && $5 >= 24000000 {print $0}' > tmp	#filter out translocation
		mv tmp $CELL_TYPE"_22_100kb".bed
fi	

cd ..
