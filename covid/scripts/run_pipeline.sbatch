#!/bin/bash

#SBATCH -A nssac_students
#SBATCH -p bii
#SBATCH -t 00:20:00
#SBATCH --output=/sfs/lustre/bahamut/scratch/jho5ze/bionets/covid/logs/pipeline_slurm-%A.out
#SBATCH --error=/sfs/lustre/bahamut/scratch/jho5ze/bionets/covid/logs/pipeline_slurm-%A.out
#SBATCH --cpus-per-task=1

set -e 

graph_out=$1
shift
uncert=$1
shift

echo "python3 /scratch/jho5ze/bionets/covid/scripts/msa_to_graph.py $graph_out $uncert $@"

python3 /scratch/jho5ze/bionets/covid/scripts/msa_to_graph.py $graph_out $uncert $@

