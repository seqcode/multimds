set -e

RES=25000
RES_KB=$(($RES/1000))

for CELLTYPE in K562 IMR90 HUVEC HMEC
do
	./get_hic_data.sh $CELLTYPE $RES
done

#python ../multimds.py hic_data/IMR90_21_${RES_KB}kb.bed hic_data/HMEC_21_${RES_KB}kb.bed
python fig6_multimds.py IMR90 HMEC
python fig7.py IMR90 $RES_KB

#python ../multimds.py hic_data/HMEC_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python fig6_multimds.py IMR90 HUVEC
python fig7.py HMEC $RES_KB

#python ../multimds.py hic_data/IMR90_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python fig6_multimds.py HMEC HUVEC
python fig7.py HUVEC $RES_KB

#python ../multimds.py hic_data/K562_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python fig6_multimds.py K562 HUVEC
python fig7.py K562 $RES_KB
