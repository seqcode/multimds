set -e

RES=100000
for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE $RES
done
./get_activity_data.sh $RES

python sup15.py
