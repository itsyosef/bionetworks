#!/bin/bash

#SBATCH -A nssac_students
#SBATCH -p bii
#SBATCH -t 0:30:00
#SBATCH --output=/sfs/lustre/bahamut/scratch/jho5ze/bionets/covid/logs/pipeline_slurm_alignment_score-%A.out
#SBATCH --error=/sfs/lustre/bahamut/scratch/jho5ze/bionets/covid/logs/pipeline_slurm_alignment_score-%A.out
#SBATCH --cpus-per-task=1

set -e 

echo "python3 /scratch/jho5ze/bionets/covid/scripts/msa_to_alignment_score.py $@"

python3 /scratch/jho5ze/bionets/covid/scripts/msa_to_alignment_score.py $@
