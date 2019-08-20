set -e

./simulate_tads.sh
	
python make_heatmap.py 1
python make_heatmap.py 2
