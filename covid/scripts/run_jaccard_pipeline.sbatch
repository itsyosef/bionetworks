#!/bin/bash

#SBATCH -A nssac_students
#SBATCH -p bii
#SBATCH -t 05:30:00
#SBATCH --output=/sfs/lustre/bahamut/scratch/jho5ze/bionets/covid/logs/pipeline_slurm_jaccard-%A.out
#SBATCH --error=/sfs/lustre/bahamut/scratch/jho5ze/bionets/covid/logs/pipeline_slurm_jaccard-%A.out
#SBATCH --cpus-per-task=1

set -e 

echo "python3 /scratch/jho5ze/bionets/covid/scripts/msa_to_jaccard.py $@"

python3 /scratch/jho5ze/bionets/covid/scripts/msa_to_jaccard.py $@

