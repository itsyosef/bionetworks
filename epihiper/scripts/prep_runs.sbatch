#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --account=nssac_covid19
#SBATCH --partition=bii
#SBATCH --time=00:30:00
#SBATCH --mem=20G

python3 prep_runs.py $SLURM_ARRAY_TASK_ID
