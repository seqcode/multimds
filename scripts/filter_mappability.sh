set -e

PREFIX=$1
RES=$2
RES_KB=$(($RES/1000))

if [ ! -e $PREFIX"_filtered".bed ]
	then
		./get_mappability.sh

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
fi
