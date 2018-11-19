set -e

./get_yeast_data.sh
python plot_relocalization.py Gal4 16 79711
