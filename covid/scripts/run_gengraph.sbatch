#!/bin/bash

#SBATCH -A nssac_students
#SBATCH -p bii
#SBATCH -t 08:00:00

set -e 
 
module load gcc openmpi mafft

echo "Moving to $1"
cd $1

echo "Running on $2 and outputting to $3"

python3 /scratch/jho5ze/bionets/covid/GenGraph/gengraphTool.py make_genome_graph --seq_file $2 --out_file_name $3

rm temp*
rm global*
rm intermediate*