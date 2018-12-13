set -e

OUTPUT_PREFIX=$1
INPUT_PREFIX=$3
RES=$2
RES_KB=$(($RES/1000))

if [ ! -d GEO ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95476/suppl/GSE95476_40kb_${INPUT_PREFIX}G1E-ER4_uninduced_iced_interchromosomal_contact_matrices.tar.gz
		tar xzvf GSE95476_40kb_${INPUT_PREFIX}G1E-ER4_uninduced_iced_interchromosomal_contact_matrices.tar.gz
		rm GSE95476_40kb_${INPUT_PREFIX}G1E-ER4_uninduced_iced_interchromosomal_contact_matrices.tar.gz
fi

for CHROM in `seq 1 19`
do
	#convert to BED
	if [ ! -e ${OUTPUT_PREFIX}-G1E_${CHROM}_40kb.bed ]
		then
			cat GEO/40kb_${INPUT_PREFIX}G1E-ER4_uninduced_iced_chr$CHROM.matrix | awk -v chrom=$CHROM '{print "chr"chrom"\t"$1"\t"$1+40000"\tchr"chrom"\t"$2"\t"$2+40000"\t"$3}' > ${OUTPUT_PREFIX}-G1E_${CHROM}_40kb.bed
	fi

	#bin BED
	if [ ! -e ${OUTPUT_PREFIX}-G1E_${CHROM}_${RES_KB}kb.bed ]
		then
			python bin_bed.py ${OUTPUT_PREFIX}-G1E_${CHROM}_40kb.bed $RES hic_data/${OUTPUT_PREFIX}-G1E_${CHROM}_${RES_KB}kb.bed
	fi
done

#convert to BED
if [ ! -e ${OUTPUT_PREFIX}-G1E_X_40kb.bed ]
	then	
		cat GEO/40kb_${INPUT_PREFIX}G1E-ER4_uninduced_iced_chrX.matrix | awk '{print "chrX\t"$1"\t"$1+40000"\tchrX\t"$2"\t"$2+40000"\t"$3}' > ${OUTPUT_PREFIX}-G1E_X_40kb.bed
fi
	
#bin BED
if [ ! -e ${OUTPUT_PREFIX}-G1E_X_${RES_KB}kb.bed ]
	then
		python bin_bed.py ${OUTPUT_PREFIX}-G1E_X_40kb.bed $RES hic_data/${OUTPUT_PREFIX}-G1E_X_${RES_KB}kb.bed
fi
