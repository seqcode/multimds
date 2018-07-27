set -e

MINIMDS_DIR=$1

./get_hic_data.sh $MINIMDS_DIR GM12878_combined
./get_hic_data.sh $MINIMDS_DIR K562

if [ ! -e peaks_filtered.bed ]
	then
		./relocalization_peaks.sh $MINIMDS_DIR
fi

if [ ! -d GM12878_combined_tadlib_input ]
	then
		python tadlib_input.py GM12878_combined
fi

if [ ! -d K562_tadlib_input ]
	then
		python tadlib_input.py K562
fi

if [ ! -e GM12878_combined_tadlib_output.txt ]
	then
		hitad -d GM12878_combined_metadata.txt -O GM12878_combined_tadlib_output.txt
fi

if [ ! -e K562_tadlib_output.txt ]
	then
		hitad -d K562_metadata.txt -O K562_tadlib_output.txt
fi

if [ ! -e GM12878_combined_K562_100kb_differential_tad_boundaries.bed ]
	then
		python differential_tad_boundaries.py
fi

#./filter_mappability.sh peaks

NUM_PEAKS=$(cat peaks_filtered.bed | wc -l)
NUM_OVERLAP=$(bedtools intersect -a GM12878_combined_K562_100kb_differential_tad_boundaries.bed -b peaks_filtered.bed | wc -l)

#negative control
if [ ! -e A_background_filtered.bed ]
	then
		python get_a_compartment.py
		bedtools subtract -a A_compartment.bed -b peaks.bed > A_background.bed
		./filter_mappability.sh A_background
fi

python tad_negative_control.py $NUM_PEAKS $NUM_OVERLAP
