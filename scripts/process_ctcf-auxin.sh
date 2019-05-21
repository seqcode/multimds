set -e

if [ ! -e GSM264494*.100000.cool ]
	then
		tar xvf GSE98671_RAW.tar
		rm GSE98671_RAW.tar
		gunzip *.gz
fi

if [ ! -e hic_data/mESC-CTCF-KO-rep1_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644947_Auxin2days-R1.100000.cool hic_data/mESC-CTCF-KO-rep1_100kb.bed
fi
if [ ! -e hic_data/mESC-CTCF-KO-rep2_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644948_Auxin2days-R2.100000.cool hic_data/mESC-CTCF-KO-rep2_100kb.bed
fi

for REP in 1 2
do
	./split_by_chrom.sh hic_data/mESC-CTCF-KO-rep${REP}_100kb.bed
done
