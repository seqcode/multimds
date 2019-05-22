set -e

for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE 100000
done

python embedding_error.py

mv GM12878_combined_K562_embedding_error.png sup1.png
