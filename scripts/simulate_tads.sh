set -e

for NUM in 1 2
do
	if [ ! -e sim${NUM}_chr21.fastq ]
		then
			sim3C --linear --create-cids --dist uniform -n 500000 -l 150 -e HindIII -m hic /data/drive1/genomes/hg19/chr21.fa sim${NUM}_chr21.fastq
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
