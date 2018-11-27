set -e

#./process_ctcf.sh

#CELLTYPES=(mESC-WT-rep1 mESC-WT-rep2 mESC-CTCF-KO-rep1 mESC-CTCF-KO-rep2)

#for i in `seq 0 $((${#CELLTYPES[@]}-1))`
#do	
#	for j in `seq 0 $(($i-1))`
#	do
#		python test_quantify_z.py ${CELLTYPES[$i]} ${CELLTYPES[$j]} 20
#	done
#done

if [ -e ctcf_design.txt ]
	then
		rm ctcf_design.txt
fi

for CELLTYPE1 in mESC-WT-rep1 mESC-WT-rep2
do
	echo $CELLTYPE1" mESC-CTCF-KO-rep1" >> ctcf_design.txt
done

python quantify_z.py 20 ctcf_design.txt 0.015
