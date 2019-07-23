set -e

RES=100000
for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE $RES
done
./get_activity_data.sh $RES

python ../multimds.py hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
python sup11.py
