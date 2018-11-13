set -e

./get_yeast_data.sh
python plot_relocalization.py Has1-Tda1 13 852000
