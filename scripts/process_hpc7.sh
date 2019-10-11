set -e

RES=$1
RES_KB=$(($RES/1000))

for REP in 1 2
do
	if [ ! -s E-MTAB-3954.HPC7_HiC_bioRep${REP}_uniques.bam ]
		then
			wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3954/E-MTAB-3954.HPC7_HiC_bioRep${REP}_uniques.bam
	fi

	for CHROM in `seq 1 19`
	do
		echo $CHROM

		#convert to BED
		if [ ! -s HPC7-rep${REP}_${CHROM}.bed ]
			then
				samtools view E-MTAB-3954.HPC7_HiC_bioRep${REP}_uniques.bam | awk -v chrom=$CHROM '$3 == chrom && $7 == "=" {print "chr"$3"\t"$4"\t"$4+1"\tchr"$3"\t"$8"\t"$8+1"\t1"}' > HPC7-rep${REP}_${CHROM}.bed
		fi
	
		#bin BED
		if [ ! -s HPC7-rep${REP}_${CHROM}_${RES_KB}kb.bed ]
			then
				python bin_bed.py HPC7-rep${REP}_${CHROM}.bed $RES HPC7-rep${REP}_${CHROM}_${RES_KB}kb.bed
		fi
	done
done
