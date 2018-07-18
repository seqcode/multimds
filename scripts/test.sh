set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

python ../multimds.py --full /data/drive1/hic_data/GM12878_combined_21_100kb.bed /data/drive1/hic_data/K562_21_100kb.bed
python ../multimds.py -N 2 /data/drive1/hic_data/GM12878_combined_21_100kb.bed /data/drive1/hic_data/K562_21_100kb.bed
python ../multimds.py -l 5 /data/drive1/hic_data/GM12878_combined_21_100kb.bed /data/drive1/hic_data/K562_21_100kb.bed
python ../multimds.py -n 1 /data/drive1/hic_data/GM12878_combined_21_100kb.bed /data/drive1/hic_data/K562_21_100kb.bed
python ../multimds.py -a 3 /data/drive1/hic_data/GM12878_combined_21_100kb.bed /data/drive1/hic_data/K562_21_100kb.bed
python ../multimds.py -m 28000000 /data/drive1/hic_data/GM12878_combined_20_100kb.bed /data/drive1/hic_data/K562_20_100kb.bed
