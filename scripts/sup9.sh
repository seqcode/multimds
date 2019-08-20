set -e

for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE 100000
done

python sup9a.py
python sup9b.py
