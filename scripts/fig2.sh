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
		cat hic_data/${CONDITION}_${STRAIN}_12_32kb.bed | awk '$2 >= 448000 && $5 >= 448000 {print $0}' > hic_data/${CONDITION}_${STRAIN}_12-downstream_32kb.bed
	done

	python plot_relocalization.py rDNA 12-downstream 448000 $STRAIN
done

for CONDITION in ctrl galactose
	do
		cat hic_data/${CONDITION}_Scer_12_32kb.bed | awk '$2 < 448000 && $5 < 448000 {print $0}' > hic_data/${CONDITION}_Scer_12-upstream_32kb.bed
done

python plot_relocalization.py Gal2 12-upstream 290212 Scer
python plot_relocalization.py Gal2 15 740086 Suva
