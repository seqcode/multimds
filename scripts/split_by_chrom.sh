set -e

BEDFILE=$1
PREFIX=${BEDFILE%.bed}

for CHROM in `seq 19`
do
	if [ ! -e ${PREFIX}_${CHROM}.bed ]
		then
			cat $BEDFILE | awk -v chrom="chr"$CHROM '$1 == chrom && $4 == chrom {print $0}' > ${PREFIX}_${CHROM}.bed
	fi
done
