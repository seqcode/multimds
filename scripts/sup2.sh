set -e

./process_ctcf-wt.sh
./process_cohesin-wt.sh

cd hic_data
python ../reproducibility.py mESC-WT-rep1_19_100kb.bed hepatocyte-WT_19_100kb.bed
cd ..
mv hic_data/GM12878_combined_K562_reproducibility.png sup2.png
