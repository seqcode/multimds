set -e

for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE 100000
done

python reproducibility.py hic_data/GM12878_combined_21_100kb.bed hic_data/K562_21_100kb.bed
mv GM12878_combined_21_100kb_K562_21_100kb_reproducibility.png fig1c.png
