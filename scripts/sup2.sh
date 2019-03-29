set -e

./get_yeast_data.sh

for STRAIN in Scer Suva
do
	python sup2.py $STRAIN 13 852000
	python sup2.py $STRAIN 2 277000
	python sup2.py $STRAIN 4 463434
	python sup2.py $STRAIN 16 79711
	
	for CONDITION in ctrl galactose
	do
		cat hic_data/${CONDITION}_${STRAIN}_12_32kb.bed | awk '$2 >= 448000 && $5 >= 448000 {print $0}' > hic_data/${CONDITION}_${STRAIN}_12-downstream_32kb.bed
	done

	python sup2.py $STRAIN 12-downstream 448000
done

for CONDITION in ctrl galactose
	do
		cat hic_data/${CONDITION}_Scer_12_32kb.bed | awk '$2 < 448000 && $5 < 448000 {print $0}' > hic_data/${CONDITION}_Scer_12-upstream_32kb.bed
done

python sup2.py Scer 12-upstream 290212
python sup2.py Suva 15 740086
