set -e

PREFIX=$1
RES=$2
RES_KB=$(($RES/1000))

#if [ ! -e $PREFIX"_filtered".bed ]
#	then
		if [ ! -e mappability.wig ]
			then
				if [ ! -e wgEncodeDukeMapabilityUniqueness35bp.bigWig ]
					then
						wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityUniqueness35bp.bigWig
				fi
				if [ ! -e bigWigToWig ]
					then
						wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
						chmod +x bigWigToWig
				fi
				./bigWigToWig wgEncodeDukeMapabilityUniqueness35bp.bigWig mappability.wig
				rm wgEncodeDukeMapabilityUniqueness35bp.bigWig
		fi

		if [ ! -e mappability_${RES_KB}kb.bed ]
			then
				python wig_to_bed.py mappability.wig $RES $(cat mappability.wig | wc -l)
		fi

		if [ ! -e mappability_${RES_KB}kb_sorted.bed ]
			then
				bedtools sort -i mappability_${RES_KB}kb.bed > mappability_${RES_KB}kb_sorted.bed
		fi

		bedtools sort -i $PREFIX.bed > $PREFIX"_sorted".bed
		NF=$(cat $PREFIX"_sorted".bed | awk '{print NF}' | head -1)
		if [ $NF -eq 3 ]
			then
				bedtools map -a $PREFIX"_sorted".bed -b mappability_${RES_KB}kb_sorted.bed -o mean -c 4 | awk '$4 >= 0.75 {print $1"\t"$2"\t"$3}' > $PREFIX"_filtered".bed	#filter out mappability < 0.75
		elif [ $NF -eq 4 ]
			then
				bedtools map -a $PREFIX"_sorted".bed -b mappability_${RES_KB}kb_sorted.bed -o mean -c 4 | awk '$5 >= 0.75 {print $1"\t"$2"\t"$3"\t"$4}' > $PREFIX"_filtered".bed	
		elif [ $NF -eq 5 ]
			then
				bedtools map -a $PREFIX"_sorted".bed -b mappability_${RES_KB}kb_sorted.bed -o mean -c 4 | awk '$6 >= 0.75 {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $PREFIX"_filtered".bed	
		fi
#fi
