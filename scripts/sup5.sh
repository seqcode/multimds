set -e

for CELLTYPE in GM12878_combined K562
do
	./get_hic_data.sh $CELLTYPE 5000
done

./get_mappability.sh 5000

if [ ! -s mappability_21_5kb.bed ]
	then
		cat mappability_5kb.bed | awk '$1 == "chr21" {print $0}' > mappability_21_5kb.bed
fi

python sup5.py
