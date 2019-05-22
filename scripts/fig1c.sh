set -e

for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE 100000
done

cd hic_data
python ../reproducibility.py GM12878_combined_21_100kb.bed K562_21_100kb.bed
cd ..
mv hic_data/GM12878_combined_K562_reproducibility.png fig1c.png
