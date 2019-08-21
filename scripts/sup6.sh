set -e

./simulate_tads.sh
	
python sup6.py 1
python sup6.py 2
