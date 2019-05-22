set -e

for CELLTYPE in K562 IMR90 HUVEC HMEC
do
	./get_hic_data.sh $CELLTYPE 100000
done

python plot_3d.py IMR90 HMEC
python plot_3d.py IMR90 HUVEC
python plot_3d.py HMEC HUVEC
python plot_3d.py K562 HUVEC
