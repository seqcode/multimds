set -e

RES=$1
RES_KB=$(($RES/1000))

mkdir -p hic_data

if [ ! -e juicer_tools.1.8.9_jcuda.0.8.jar ]
	then
		curl http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.8.9_jcuda.0.8.jar -o juicer_tools.1.8.9_jcuda.0.8.jar
fi

for f in *.hic
do
	for CHROM in `seq 19`
	do
		if [ ! -e hic_data/${f%.*}_${CHROM}_${RES_KB}kb.bed ]
			then 
				echo $CHROM
				OUT=${f%.*}_${CHROM}_${RES_KB}kb.tsv
				java -jar juicer_tools.1.8.9_jcuda.0.8.jar dump observed KR $f $CHROM $CHROM BP $RES $OUT
				if [ -s $OUT ]
					then
						cat $OUT | awk -v chr=$CHROM -v res=$RES '$3 != "NaN" {print "chr"chr"\t"$1"\t"$1+res"\tchr"chr"\t"$2"\t"$2+res"\t"$3}' > hic_data/${f%.*}_${CHROM}_${RES_KB}kb.bed
				fi
		fi
	done
done

rm *.tsv
