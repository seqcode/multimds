set -e

if [ ! -s GSE93431_NIPBL.100kb.cool.HDF5 ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93431/suppl/GSE93431_NIPBL.100kb.cool.HDF5.gz
		gunzip GSE93431_NIPBL.100kb.cool.HDF5.gz
fi

if [ ! -s hic_data/hepatocyte-cohesin-KO_100kb.bed ]
	then
		python bed_from_hdf5.py GSE93431_NIPBL.100kb.cool.HDF5 hic_data/hepatocyte-cohesin-KO_100kb.bed
fi

./split_by_chrom.sh hic_data/hepatocyte-cohesin-KO_100kb.bed
