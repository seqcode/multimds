set -e

if [ ! -e chr21.fa ]
	then
		curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz -o chr21.fa.gz
		gunzip chr21.fa.gz
fi

for NUM in 1 2
do
	if [ -e profile.tsv ]
		then
			rm profile.tsv
	fi

	if [ ! -e sim${NUM}_chr21.fastq ]
		then
			sim3C --linear --create-cids --dist uniform -n 1000000 -l 150 -e HindIII -m hic --depths depths$NUM.tsv --subseqs subseqs$NUM.tsv chr21.fa sim${NUM}_chr21.fastq
	fi

	if [ ! -e sim${NUM}_chr21.bed ]
		then
			python make_bed.py $NUM
	fi

	if [ ! -e sim${NUM}_chr21_100kb.bed ]
		then
			python bin_bed.py sim${NUM}_chr21.bed 100000 sim${NUM}_chr21_100kb.bed
	fi
done
