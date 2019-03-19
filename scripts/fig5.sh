set -e

#cohesin
./process_cohesin-wt.sh
./process_cohesin-ko.sh
echo "hepatocyte-cohesin-KO hepatocyte-WT" > cohesin_design.txt
python quantify_z.py 20 cohesin_design.txt 0.04 0.75 cohesin depletion

#Brd2

RES=100000
./process_g1e.sh WT $RES
./process_g1e.sh KO-rep1 $RES D_BRD2KO_
./process_g1e.sh KO-rep2 $RES C8_BRD2KO_

echo "KO-rep1-G1E KO-rep2-G1E" > brd2_rep_design.txt
echo "KO-rep1-G1E WT-G1E" > brd2_design.txt
echo "KO-rep2-G1E WT-G1E" >> brd2_design.txt

python quantify_z.py 20 brd2_design.txt 0.02 0.55 Brd2 depletion
python quantify_z.py 20 brd2_rep_design.txt 0.035

#CTCF
./process_ctcf-wt.sh
./process_ctcf-ko.sh

if [ -e ctcf_design.txt ]
	then
		rm ctcf_design.txt
fi

for CELLTYPE1 in mESC-WT-rep1 mESC-WT-rep2
do
	echo $CELLTYPE1" mESC-CTCF-KO-rep1" >> ctcf_design.txt
done

python quantify_z.py 20 ctcf_design.txt 0.015 0.75 CTCF depletion
