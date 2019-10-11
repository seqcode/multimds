set -e

if [ ! -s GSM264494*.100000.cool ]
	then
		tar xvf GSE98671_RAW.tar
		rm GSE98671_RAW.tar
		gunzip *.gz
fi

if [ ! -s hic_data/mESC-WT-rep1_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644945_Untreated-R1.100000.cool hic_data/mESC-WT-rep1_100kb.bed
fi
if [ ! -s hic_data/mESC-WT-rep2_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644946_Untreated-R2.100000.cool hic_data/mESC-WT-rep2_100kb.bed
fi

for REP in 1 2
do
	./split_by_chrom.sh hic_data/mESC-WT-rep${REP}_100kb
done
