set -e

PREFIX=$1

if [ ! -e $PREFIX"_filtered".bed ]
	then
		if [ ! -e mappability_sorted.bed ]
			then
				wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityUniqueness35bp.bigWig
				wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
				chmod +x bigWigToWig
				./bigWigToWig wgEncodeDukeMapabilityUniqueness35bp.bigWig mappability.wig
				rm wgEncodeDukeMapabilityUniqueness35bp.bigWig
				python wig_to_bed.py mappability.wig 100000 $(cat mappability.wig | wc -l)
				bedtools sort -i mappability.bed > mappability_sorted.bed
		fi
		bedtools sort -i $PREFIX.bed > $PREFIX"_sorted".bed
		NF=$(cat $PREFIX"_sorted".bed | awk '{print NF}' | head -1)
		if [ $NF -eq 3 ]
			then
				bedtools map -a $PREFIX"_sorted".bed -b mappability_sorted.bed -o mean -c 4 | awk '$4 >= 0.75 {print $1"\t"$2"\t"$3}' > $PREFIX"_filtered".bed	#filter out mappability < 0.75
		elif [ $NF -eq 4 ]
			then
				bedtools map -a $PREFIX"_sorted".bed -b mappability_sorted.bed -o mean -c 4 | awk '$5 >= 0.75 {print $1"\t"$2"\t"$3"\t"$4}' > $PREFIX"_filtered".bed	
		fi
fi
