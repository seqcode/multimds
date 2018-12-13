set -e

#./get_lymphoblastoid.sh

if [ -e lymphoblastoid_independent_design.txt ]
	then
		rm lymphoblastoid_independent_design.txt
fi

CELLTYPES=(GM19238 GM19239 GM19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733)

for i in `seq 0 $((${#CELLTYPES[@]}-1))`
do
	CELLTYPE1=${CELLTYPES[$i]}
	for j in `seq 0 $(($i-1))`
	do
		CELLTYPE2=${CELLTYPES[$j]}
		echo $CELLTYPE1" "$CELLTYPE2 >> lymphoblastoid_independent_design.txt
	done

done

python quantify_z_independent.py 20 lymphoblastoid_independent_design.txt 0.03
