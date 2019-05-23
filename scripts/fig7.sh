set -e

RES=25000
RES_KB=$(($RES/1000))

for CELLTYPE in K562 IMR90 HUVEC HMEC
do
	./get_hic_data.sh $CELLTYPE $RES
done

python ../multimds.py hic_data/IMR90_21_${RES_KB}kb.bed hic_data/HMEC_21_${RES_KB}kb.bed
python plot_3d.py IMR90 HMEC $RES_KB

python ../multimds.py hic_data/IMR90_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python plot_3d.py IMR90 HUVEC $RES_KB

python ../multimds.py hic_data/HMEC_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python plot_3d.py HMEC HUVEC $RES_KB

python ../multimds.py hic_data/K562_21_${RES_KB}kb.bed hic_data/HUVEC_21_${RES_KB}kb.bed
python plot_3d.py K562 HUVEC $RES_KB
