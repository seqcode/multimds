set -e

./get_yeast_data.sh

for STRAIN in Scer Suva
do
	python ../minimds.py -w 0 -o independent_ctrl_${STRAIN}_13_32kb_structure.tsv hic_data/ctrl_${STRAIN}_13_32kb.bed
	python ../minimds.py -w 0 -o independent_galactose_${STRAIN}_13_32kb_structure.tsv hic_data/galactose_${STRAIN}_13_32kb.bed
	python plot_relocalization.py Has1-Tda1 13 852000 $STRAIN independent_

	python ../minimds.py -w 0 -o independent_ctrl_${STRAIN}_2_32kb_structure.tsv hic_data/ctrl_${STRAIN}_2_32kb.bed
	python ../minimds.py -w 0 -o independent_galactose_${STRAIN}_2_32kb_structure.tsv hic_data/galactose_${STRAIN}_2_32kb.bed
	python plot_relocalization.py Gal1-7-10 2 277000 $STRAIN independent_

	python ../minimds.py -w 0 -o independent_ctrl_${STRAIN}_4_32kb_structure.tsv hic_data/ctrl_${STRAIN}_4_32kb.bed
	python ../minimds.py -w 0 -o independent_galactose_${STRAIN}_4_32kb_structure.tsv hic_data/galactose_${STRAIN}_4_32kb.bed
	python plot_relocalization.py Gal3 4 463434 $STRAIN independent_

	python ../minimds.py -w 0 -o independent_ctrl_${STRAIN}_16_32kb_structure.tsv hic_data/ctrl_${STRAIN}_16_32kb.bed
	python ../minimds.py -w 0 -o independent_galactose_${STRAIN}_16_32kb_structure.tsv hic_data/galactose_${STRAIN}_16_32kb.bed
	python plot_relocalization.py Gal4 16 79711 $STRAIN independent_
	
	for CONDITION in ctrl galactose
	do
		cat ${CONDITION}_${STRAIN}_12_32kb.bed | awk '$2 >= 448000 && $5 >= 448000 {print $0}' > ${CONDITION}_${STRAIN}_12-downstream_32kb.bed
	done
	
	python ../minimds.py -w 0 -o independent_ctrl_${STRAIN}_12-downstream_32kb_structure.tsv hic_data/ctrl_${STRAIN}_12-downstream_32kb.bed
	python ../minimds.py -w 0 -o independent_galactose_${STRAIN}_12-downstream_32kb_structure.tsv hic_data/galactose_${STRAIN}_12-downstream_32kb.bed
	python plot_relocalization.py rDNA 12-downstream 448000 $STRAIN independent_
done

for CONDITION in ctrl galactose
do
	cat ${CONDITION}_Scer_12_32kb.bed | awk '$2 < 448000 && $5 < 448000 {print $0}' > ${CONDITION}_Scer_12-upstream_32kb.bed
done

python ../minimds.py -w 0 -o independent_ctrl_Scer_12-upstream_32kb_structure.tsv hic_data/ctrl_Scer_12-upstream_32kb.bed
python ../minimds.py -w 0 -o independent_galactose_Scer_12-upstream_32kb_structure.tsv hic_data/galactose_Scer_12-upstream_32kb.bed
python plot_relocalization.py Gal2 12-upstream 290212 Scer independent_

python ../minimds.py -w 0 -o independent_ctrl_Suva_15_32kb_structure.tsv hic_data/ctrl_Suva_15_32kb.bed
python ../minimds.py -w 0 -o independent_galactose_Suva_15_32kb_structure.tsv hic_data/galactose_Suva_15_32kb.bed
python plot_relocalization.py Gal2 15 740086 Suva independent_
