set -e

#./process_cohesin.sh
echo "hepatocyte-cohesin-KO hepatocyte-WT" > cohesin_independent_design.txt
python quantify_z_independent.py 20 cohesin_independent_design.txt 0.04
