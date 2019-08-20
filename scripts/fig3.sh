set -e

./get_yeast_data.sh

for STRAIN in Scer Suva
do
	for CONDITION in ctrl galactose
	do
		cat hic_data/${CONDITION}_${STRAIN}_12_32kb.bed | awk '$2 >= 448000 && $5 >= 448000 {print $0}' > hic_data/${CONDITION}_${STRAIN}_12-downstream_32kb.bed
	done
	
	for CHROM in 13 2 4 16 12-downstream
	do
		python fig3_multimds.py $CHROM $STRAIN
	done

	#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_${STRAIN}_13_32kb.bed hic_data/galactose_${STRAIN}_13_32kb.bed
	python fig3.py Has1-Tda1 13 852000 $STRAIN

	#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_${STRAIN}_2_32kb.bed hic_data/galactose_${STRAIN}_2_32kb.bed
	python fig3.py Gal1-7-10 2 277000 $STRAIN

	#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_${STRAIN}_4_32kb.bed hic_data/galactose_${STRAIN}_4_32kb.bed 
	python fig3.py Gal3 4 463434 $STRAIN

	#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_${STRAIN}_16_32kb.bed hic_data/galactose_${STRAIN}_16_32kb.bed
	python fig3.py Gal4 16 79711 $STRAIN
	
	#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_${STRAIN}_12-downstream_32kb.bed hic_data/galactose_${STRAIN}_12-downstream_32kb.bed
	python fig3.py rDNA 12-downstream 448000 $STRAIN
done

for CONDITION in ctrl galactose
do
	cat hic_data/${CONDITION}_Scer_12_32kb.bed | awk '$2 < 448000 && $5 < 448000 {print $0}' > hic_data/${CONDITION}_Scer_12-upstream_32kb.bed
done

python fig3_multimds.py 12-upstream Scer
#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_Scer_12-upstream_32kb.bed hic_data/galactose_Scer_12-upstream_32kb.bed
python fig3.py Gal2 12-upstream 290212 Scer

python fig3_multimds.py 15 Suva
#python ../multimds.py -P 0.1 -w 0 hic_data/ctrl_Suva_15_32kb.bed hic_data/galactose_Suva_15_32kb.bed
python fig3.py Gal2 15 740086 Suva
