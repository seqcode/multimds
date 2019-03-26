set -e

if [ -e output ]
	then
		rm output
fi

for CHROM in 21 20 19 18 17 16 15 14 13 12 11 10 8 7 6 5 4 3 2 1
do
	python get_chrom_size.py $CHROM >> output
done

cat output | awk '$1 == "size" {print $2}' > chrom_sizes.txt
rm output
