set -e

for CELLTYPE in GM12878_primary GM12878_replicate
do
	./get_hic_data.sh $CELLTYPE 100000
done

echo "GM12878_primary GM12878_replicate" > encode_rep_design.txt

python quantify_z.py 23 encode_rep_design.txt 0.025 0.5 GM12878 reps
