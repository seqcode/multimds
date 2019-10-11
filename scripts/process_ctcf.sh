set -e

if [ ! -s GSM264494*.100000.cool ]
	then
		tar xvf GSE98671_RAW.tar
		rm GSE98671_RAW.tar
		gunzip *.gz
fi

if [ ! -s mESC-WT-rep1_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644945_Untreated-R1.100000.cool hic_data/mESC-WT-rep1_100kb.bed
fi
if [ ! -s mESC-WT-rep2_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644946_Untreated-R2.100000.cool hic_data/mESC-WT-rep2_100kb.bed
fi
if [ ! -s mESC-CTCF-KO-rep1_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644947_Auxin2days-R1.100000.cool hic_data/mESC-CTCF-KO-rep1_100kb.bed
fi
if [ ! -s mESC-CTCF-KO-rep2_100kb.bed ]
	then
		python bed_from_hdf5.py GSM2644948_Auxin2days-R2.100000.cool hic_data/mESC-CTCF-KO-rep2_100kb.bed
fi

for PREFIX in mESC-WT-rep1 mESC-WT-rep2 mESC-CTCF-KO-rep1 mESC-CTCF-KO-rep2
do
	./split_by_chrom.sh hic_data/$PREFIX
done
