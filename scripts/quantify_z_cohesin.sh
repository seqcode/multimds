set -e

#./process_cohesin.sh
echo "hepatocyte-cohesin-KO hepatocyte-WT" > cohesin_design.txt
python quantify_z.py 20 cohesin_design.txt 0.04
