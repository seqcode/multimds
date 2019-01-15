set -e

for GENE in Gal1-7-10 Gal2 Has1-Tda1 Gal3 Gal4 Hxt1
do
	python pileup.py $GENE
done
