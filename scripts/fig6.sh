set -e

RES=100000
RES_KB=$(($RES/1000))

for CELLTYPE in K562 IMR90 HUVEC HMEC
do
	./get_hic_data.sh $CELLTYPE $RES
done

python ../multimds.py hic_data/IMR90_21_${RES_KB}kb.bed hic_data/HMEC_21_${RES_KB}kb.bed
python fig6_multimds.py IMR90 HMEC
python fig6.py IMR90 HMEC $RES_KB

python ../multimds.py hic_data/IMR90_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python fig6_multimds.py IMR90 HUVEC
python fig6.py IMR90 HUVEC $RES_KB

python ../multimds.py hic_data/HMEC_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python fig6_multimds.py HMEC HUVEC
python fig6.py HMEC HUVEC $RES_KB

python ../multimds.py hic_data/K562_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python fig6_multimds.py K562 HUVEC
python fig6.py K562 HUVEC $RES_KB
