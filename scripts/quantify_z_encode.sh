set -e

for CELLTYPE in GM12878_primary GM12878_replicate K562 KBM7 IMR90 HUVEC HMEC NHEK
do
	./get_hic_data.sh $CELLTYPE
done

if [ -e encode_design.txt ]
	then
		rm encode_design.txt
fi

for CELLTYPE1 in GM12878_primary GM12878_replicate
do
	for CELLTYPE2 in K562 IMR90 HUVEC HMEC NHEK
	do
		echo $CELLTYPE1" "$CELLTYPE2 >> encode_design.txt
	done
done

CELLTYPES=(K562 IMR90 HUVEC HMEC NHEK)

for i in `seq 0 $((${#CELLTYPES[@]}-1))`
do	
	for j in `seq 0 $(($i-1))`
	do
		echo ${CELLTYPES[$i]}" "${CELLTYPES[$j]} >> encode_design.txt
	done
done

echo "GM12878_primary GM12878_replicate" > encode_rep_design.txt

python quantify_z.py 23 encode_design.txt 0.025
python quantify_z.py 23 encode_rep_design.txt 0.025
