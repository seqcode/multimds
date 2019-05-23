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

#align
for NUM in `seq 3 6`
do
	SRR=545748$NUM
 
	if [ ! -e SRR${SRR}.sam ]
		then

			#download
			curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR545/SRR$SRR/SRR$SRR.sra -o SRR$SRR.sra

			#fastq-dump
			fastq-dump --split-files SRR$SRR.sra
			rm SRR$SRR.sra

		#index
		for READ in 1 2
		do
			bwa aln /data/drive1/genomes/Sc_Su/Sc_Su.fa SRR${SRR}_$READ.fastq > SRR${SRR}_$READ.sai
		done

		#map
		bwa sampe /data/drive1/genomes/Sc_Su/Sc_Su.fa SRR${SRR}_1.sai SRR${SRR}_2.sai SRR${SRR}_1.fastq SRR${SRR}_2.fastq > SRR${SRR}.sam
	fi
done

#call peaks independently
if [ ! -e galactose_Nup60_peaks.broadPeak ]
	then
		macs2 callpeak --broad -t SRR5457483.sam -c SRR5457484.sam -n galactose_Nup60 --nomodel
fi
if [ ! -e glucose_Nup60_peaks.broadPeak ]
	then
		macs2 callpeak --broad -t SRR5457485.sam -c SRR5457486.sam -n glucose_Nup60 --nomodel
fi

#merge
cp galactose_Nup60_peaks.broadPeak all_Nup60.bed
cat glucose_Nup60_peaks.broadPeak >> all_Nup60.bed

bedtools sort -i all_Nup60.bed > all_Nup60_sorted.bed
bedtools merge -i all_Nup60_sorted.bed > all_Nup60_merged.bed

#convert to bam
for NUM in 3 5
do
	if [ ! -e SRR545748$NUM.bam ]
		then
			samtools view -b SRR545748$NUM.sam > SRR545748$NUM.bam
	fi
done

#get tag counts
bedtools coverage -a all_Nup60_merged.bed -b SRR5457483.bam -counts > galactose_coverage.bed
bedtools coverage -a all_Nup60_merged.bed -b SRR5457485.bam -counts > glucose_coverage.bed

echo "Symbol	ctrl	galactose" > nup60_counts.tsv
paste glucose_coverage.bed galactose_coverage.bed | cut -d "	" -f 1,2,3,4,8 | awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' >> nup60_counts.tsv

#get differential peaks
Rscript run_edger.R
python process_edger_results.py

for GENE in Gal1-7-10 Gal2 Has1-Tda1 Gal3 Gal4 Hxt1
do
	python pileup.py $GENE
done
