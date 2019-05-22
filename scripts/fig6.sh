set -e

#ENCODE
for CELLTYPE in K562 IMR90 HUVEC HMEC
do
	./get_hic_data.sh $CELLTYPE 100000
done

python plot_relocalization_compartments.py IMR90 HMEC
python plot_relocalization_compartments.py IMR90 HUVEC
python plot_relocalization_compartments.py HMEC HUVEC
python plot_relocalization_compartments.py K562 HUVEC
