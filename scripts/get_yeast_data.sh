set -e

if [ ! -e ctrl_32kb_matrix.txt.gz ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2355nnn/GSM2355770/suppl/GSM2355770_ILY456_exponential_Sau3AI.ScSu.32000.matrix.txt.gz -o ctrl_32kb_matrix.txt.gz
		gunzip ctrl_32kb_matrix.txt.gz
fi

if [ ! -e galactose_32kb_matrix.txt.gz ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2355nnn/GSM2355772/suppl/GSM2355772_ILY456_galactose_Sau3AI.ScSu.32000.matrix.txt.gz -o galactose_32kb_matrix.txt.gz
		gunzip galactose_32kb_matrix.txt.gz
fi

if [ ! -e GSE88952_Sc_Su.32000.bed ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn/GSE88952/suppl/GSE88952_Sc_Su.32000.bed.gz -o GSE88952_Sc_Su.32000.bed.gz
		gunzip GSE88952_Sc_Su.32000.bed.gz
fi

python convert_to_bed.py

for PREFIX in ctrl galactose
do
	for STRAIN in Scer Suva
	do
		for CHROM_NUM in `seq 16`
		do
			CHROM_NAME=${STRAIN}_$CHROM_NUM
			if [ ! -e ${PREFIX}_${CHROM_NAME}_32kb.bed ]
				then
					echo $CHROM_NAME
					cat ${PREFIX}_32kb.bed | awk -v chrom_name=$CHROM_NAME '$1 == chrom_name && $4 == chrom_name && $7 != 0 && $7 != "nan" {print $0}' > ${PREFIX}_${CHROM_NAME}_32kb.bed
			fi
		done
	done
done
