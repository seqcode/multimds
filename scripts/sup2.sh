set -e

./process_ctcf-wt.sh
./process_cohesin-wt.sh

python reproducibility.py hic_data/mESC-WT-rep1_19_100kb.bed hic_data/hepatocyte-WT_19_100kb.bed
mv mESC-WT-rep1_19_100kb_hepatocyte-WT_19_100kb_reproducibility.png sup2.png
