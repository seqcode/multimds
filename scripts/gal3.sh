set -e

./get_yeast_data.sh
python plot_relocalization.py Gal3 4 463434
