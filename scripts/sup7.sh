set -e

#TODO: download hg19 chr21
./simulate_tads.sh
python sup7.py
