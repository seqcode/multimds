set -e

MINIMDS_DIR=$1

mkdir -p hic_data

cd hic_data

if [ ! -e GM12878_combined_*_100kb.bed ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz

		for CHROM in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22	#skip chr9
		do
			tar xzf GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz GM12878_combined/100kb_resolution_intrachromosomal/chr$CHROM
			python $MINIMDS_DIR/scripts/normalize.py GM12878_combined 100000 $CHROM
		done
fi	


if [ ! -e K562_*_100kb.bed ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_intrachromosomal_contact_matrices.tar.gz

		for CHROM in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22	#skip chr9
		do
			tar xzf GSE63525_K562_intrachromosomal_contact_matrices.tar.gz K562/100kb_resolution_intrachromosomal/chr$CHROM
			python $MINIMDS_DIR/scripts/normalize.py K562 100000 $CHROM
		done
fi

cd ..
