set -e

PREFIX=$1

cat $PREFIX.bed | awk -v start=$2 -v end=$3 '($2 < start || $2 > end) && ($5 < start || $5 > end) {print $0}' > ${PREFIX}_filtered.bed
