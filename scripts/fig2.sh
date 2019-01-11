set -e

./get_yeast_data.sh

for STRAIN in Scer Suva
do
	python plot_relocalization.py Has1-Tda1 13 852000 $STRAIN
	python plot_relocalization.py Gal1-7-10 2 277000 $STRAIN
	python plot_relocalization.py Gal3 4 463434 $STRAIN
	python plot_relocalization.py Gal4 16 79711 $STRAIN
	
	for CONDITION in ctrl galactose
	do
		cat ${CONDITION}_${STRAIN}_12_32kb.bed | awk '$2 >= 448000 && $5 >= 448000 {print $0}' > ${CONDITION}_${STRAIN}_12-downstream_32kb.bed
	done

	python plot_relocalization.py 12-downstream rDNA 448000 $STRAIN
done

for CONDITION in ctrl galactose
	do
		cat ${CONDITION}_Scer_12_32kb.bed | awk '$2 < 448000 && $5 < 448000 {print $0}' > ${CONDITION}_Scer_12-upstream_32kb.bed
done

python plot_relocalization.py 12-upstream Gal2 290212 Scer
python plot_relocalization.py 15 Gal2 740086 Suva
