set -e

for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE 100000
done

python reproducibility.py /data/drive1/hic_data/GM12878_combined_21_100kb.bed /data/drive1/hic_data/K562_21_100kb.bed
mv GM12878_combined_21_100kb_K562_21_100kb_reproducibility.svg fig1c.svg
