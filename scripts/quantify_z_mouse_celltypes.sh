set -e

RES=100000

./process_hpc7.sh $RES
#./process_g1e.sh WT
#TODO: process other cell types

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

python quantify_z.py 20 mouse_celltype_design.txt 0.035
#python quantify_z.py 20 mouse_celltype_rep_design.txt 0.02
