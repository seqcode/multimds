set -e

./get_yeast_data.sh
python plot_relocalization.py Gal1-7-10 2 277000
