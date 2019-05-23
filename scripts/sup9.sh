set -e

if [ ! -e GSM2583734_galactose_IP.bedgraph ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583734/suppl/GSM2583734_galactose_IP.bedgraph.gz -o GSM2583734_galactose_IP.bedgraph.gz
		gunzip GSM2583734_galactose_IP.bedgraph.gz
fi

if [ ! -e GSM2583735_galactose_input.bedgraph ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583735/suppl/GSM2583735_galactose_input.bedgraph.gz -o GSM2583735_galactose_input.bedgraph.gz
		gunzip GSM2583735_galactose_input.bedgraph.gz
fi

if [ ! -e GSM2583736_glucose_IP.bedgraph ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583736/suppl/GSM2583736_glucose_IP.bedgraph.gz -o GSM2583736_glucose_IP.bedgraph.gz
		gunzip GSM2583736_glucose_IP.bedgraph.gz
fi

if [ ! -e GSM2583737_glucose_input.bedgraph ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583737/suppl/GSM2583737_glucose_input.bedgraph.gz -o GSM2583737_glucose_input.bedgraph.gz
		gunzip GSM2583737_glucose_input.bedgraph.gz
fi

python sup9a.py

if [ ! -e GSM2583738_asy_r1.counts.txt ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583738/suppl/GSM2583738_asy_r1.counts.txt.gz -o GSM2583738_asy_r1.counts.txt.gz
		gunzip GSM2583738_asy_r1.counts.txt.gz
fi

if [ ! -e GSM2583739_asy_r2.counts.txt ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583739/suppl/GSM2583739_asy_r2.counts.txt.gz -o GSM2583739_asy_r2.counts.txt.gz
		gunzip GSM2583739_asy_r2.counts.txt.gz
fi

if [ ! -e GSM2583740_asy_r3.counts.txt ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583740/suppl/GSM2583740_asy_r3.counts.txt.gz	-o GSM2583740_asy_r3.counts.txt.gz
		gunzip GSM2583740_asy_r3.counts.txt.gz
fi

if [ ! -e GSM2583741_gal_r1.counts.txt ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583741/suppl/GSM2583741_gal_r1.counts.txt.gz -o GSM2583741_gal_r1.counts.txt.gz
		gunzip GSM2583741_gal_r1.counts.txt.gz
fi

if [ ! -e GSM2583742_gal_r2.counts.txt ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583742/suppl/GSM2583742_gal_r2.counts.txt.gz -o GSM2583742_gal_r2.counts.txt.gz
		gunzip GSM2583742_gal_r2.counts.txt.gz
fi

if [ ! -e GSM2583743_gal_r3.counts.txt ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583743/suppl/GSM2583743_gal_r3.counts.txt.gz -o GSM2583743_gal_r3.counts.txt.gz
		gunzip GSM2583743_gal_r3.counts.txt.gz
fi

paste GSM2583738_asy_r1.counts.txt  GSM2583739_asy_r2.counts.txt  GSM2583740_asy_r3.counts.txt  GSM2583741_gal_r1.counts.txt  GSM2583742_gal_r2.counts.txt  GSM2583743_gal_r3.counts.txt | cut -d "	" -f 1,2,4,6,8,10,12 >> rnaseq_counts.tsv

python sup9b.py
