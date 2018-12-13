set -e

#RES=100000
#./process_g1e.sh WT $RES
#./process_g1e.sh KO-rep1 $RES D_BRD2KO_
#./process_g1e.sh KO-rep2 $RES C8_BRD2KO_

echo "KO-rep1-G1E KO-rep2-G1E" > brd2_rep_independent_design.txt
echo "KO-rep1-G1E WT-G1E" > brd2_independent_design.txt
echo "KO-rep2-G1E WT-G1E" >> brd2_independent_design.txt

python quantify_z.py 20 brd2_independent_design.txt 0.02
python quantify_z.py 20 brd2_rep_independent_design.txt 0.035
