set -e

./get_hic_data.sh GM12878_combined
./get_hic_data.sh K562

for METHOD in MultiMDS Independent_MDS
do
	OUTPUT=$METHOD"_times".txt

	if [ -e $OUTPUT ]
		then
			rm $OUTPUT
	fi

	for CHROM in 22 21 20 19 18 17 16 15 14 13 12 11 10 8 7 6 5 4 3 2 1
	do
		echo $CHROM
		/usr/bin/time -o $OUTPUT -a -f %e python test_$METHOD.py $CHROM | awk '$1 == "rmsd" {print $2}' >> $OUTPUT
	done

done

if [ ! -e chrom_sizes.txt ]
	then
		./get_chrom_sizes.sh
fi

python time_methods.py
