set -e

#ENCODE
for CELLTYPE in GM12878_primary GM12878_replicate K562 KBM7 IMR90 HUVEC HMEC NHEK
do
	./get_hic_data.sh $CELLTYPE 100000
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

#mouse cell types
RES=100000

./process_hpc7.sh $RES
./process_g1e.sh WT $RES
./process_g1e.sh KO $RES
./process_ctcf-wt.sh
./process_cohesin-wt.sh
./process_cohesin-ko.sh

echo "hepatocyte-cohesin-KO hepatocyte-WT" > cohesin_design.txt

if [ -e mouse_celltype_design.txt ]
	then
		rm mouse_celltype_design.txt
fi

CELLTYPES=(mESC-WT-rep1 mESC-WT-rep2 HPC7-rep1 HPC7-rep2 WT-G1E hepatocyte-WT)

for CELLTYPE1 in mESC-WT-rep1 mESC-WT-rep2
do
	for CELLTYPE2 in HPC7-rep1 HPC7-rep2 WT-G1E hepatocyte-WT
	do
		echo $CELLTYPE1" "$CELLTYPE2 >> mouse_celltype_design.txt
	done
done

for CELLTYPE1 in HPC7-rep1 HPC7-rep2
do
	for CELLTYPE2 in mESC-WT-rep1 mESC-WT-rep2 WT-G1E hepatocyte-WT
	do
		echo $CELLTYPE1" "$CELLTYPE2 >> mouse_celltype_design.txt
	done
done

echo "WT-G1E hepatocyte-WT" >> mouse_celltype_design.txt

echo "mESC-WT-rep1 mESC-WT-rep2" > mouse_celltype_rep_design.txt
echo "HPC7-rep1 HPC7-rep2" >> mouse_celltype_rep_design.txt


#LCL
./process_lymphoblastoid.sh

if [ -e lymphoblastoid_design.txt ]
	then
		rm lymphoblastoid_design.txt
fi

CELLTYPES=(GM19238 GM19239 GM19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733)

for i in `seq 0 $((${#CELLTYPES[@]}-1))`
do
	CELLTYPE1=${CELLTYPES[$i]}
	for j in `seq 0 $(($i-1))`
	do
		CELLTYPE2=${CELLTYPES[$j]}
		echo $CELLTYPE1" "$CELLTYPE2 >> lymphoblastoid_design.txt
	done

done

python comp_correlation.py
