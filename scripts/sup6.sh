set -e

curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583734/suppl/GSM2583734_galactose_IP.bedgraph.gz -o GSM2583734_galactose_IP.bedgraph.gz
curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583735/suppl/GSM2583735_galactose_input.bedgraph.gz -o GSM2583735_galactose_input.bedgraph.gz
curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583736/suppl/GSM2583736_glucose_IP.bedgraph.gz -o GSM2583736_glucose_IP.bedgraph.gz
curl ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2583nnn/GSM2583737/suppl/GSM2583737_glucose_input.bedgraph.gz -o GSM2583737_glucose_input.bedgraph.gz

gunzip *.bedgraph.gz

#align
for NUM in `seq 3 6`
do
	SRR=SRR545748$NUM
 
	#download
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR545/SRR$SRR/SRR$SRR.sra

	#fastq-dump
	fastq-dump --split-files SRR$SRR.sra
	rm SRR$SRR.sra

	#aln
	for READ in 1 2
	do
		bwa aln /data/drive1/genomes/Sc_Su/Sc_Su.fa SRR${SRR}_$READ.fastq > SRR${SRR}_$READ.sai
	done

	#sampe
	bwa sampe /data/drive1/genomes/Sc_Su/Sc_Su.fa SRR${SRR}_1.sai SRR${SRR}_2.sai SRR${SRR}_1.fastq SRR${SRR}_2.fastq > SRR${SRR}.sam
done

#call peaks independently
macs2 callpeak --broad -t SRR5457483.sam -c SRR5457484.sam -n galactose_Nup60
macs2 callpeak --broad -t SRR5457485.sam -c SRR5457486.sam -n glucose_Nup60

#merge
cp galactose_Nup60.broadPeak all_Nup60.bed
cat glucose_Nup60.broadPeak >> all_Nup60.bed

bedtools merge -i all_Nup60.bed > all_Nup60_merged.bed

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
