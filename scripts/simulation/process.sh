set -e

for NUM in 1 2
do
	#TODO: run sim3c
	python make_bed.py $NUM
	python bin_bed.py sim${NUM}_chr21.bed 100000 sim${NUM}_chr21_100kb.bed
done
	
python make_heatmap.py 1
python make_heatmap.py 2
