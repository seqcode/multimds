set -e

./get_yeast_data.sh

for CONDITION in ctrl galactose
do
	cat ${CONDITION}_Scer_12_32kb.bed | awk '$2 < 448000 && $5 < 448000 {print $0}' > ${CONDITION}_Scer_12-upstream_32kb.bed
	cat ${CONDITION}_Scer_12_32kb.bed | awk '$2 >= 448000 && $5 >= 448000 {print $0}' > ${CONDITION}_Scer_12-downstream_32kb.bed
done

python sup4.py
