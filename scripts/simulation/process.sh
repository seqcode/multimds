set -e

for NUM in 1 2
do
	sim3C --dist uniform -n 500000 -l 150 -e HindIII -m hic /data/drive1/genomes/hg19/chr21.fa sim${NUM}_chr21.fastq
	python make_bed.py $NUM
	python bin_bed.py sim${NUM}_chr21.bed 100000 sim${NUM}_chr21_100kb.bed
done
	
python make_heatmap.py 1
python make_heatmap.py 2
