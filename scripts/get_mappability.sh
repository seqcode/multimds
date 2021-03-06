set -e

RES=$1
RES_KB=$(($RES/1000))

if [ ! -s mappability.wig ]
	then
		if [ ! -s wgEncodeDukeMapabilityUniqueness35bp.bigWig ]
			then
				curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityUniqueness35bp.bigWig -o mappability.bigWig
		fi
		./bigWigToWig mappability.bigWig mappability.wig
		rm mappability.bigWig
fi

if [ ! -s mappability_${RES_KB}kb.bed ]
	then
		python wig_to_bed.py mappability.wig $RES $(cat mappability.wig | wc -l)
fi
