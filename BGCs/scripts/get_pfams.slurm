#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH -A nssac_students
#SBATCH -p bii
#SBATCH -c 4
#SBATCH --output=get_pfams_%a.out

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    SLURM_ARRAY_TASK_ID=0
fi

gids=($(cat /sfs/lustre/bahamut/scratch/jho5ze/bionets/BGCs/misc_data/all_anti_gids_more_than_5_bgcs.txt))

gid=${gids[$SLURM_ARRAY_TASK_ID]}

profile=/sfs/lustre/bahamut/scratch/jho5ze/bionets/BGCs/misc_data/Pfam-A.hmm

#wget ftp://ftp.patricbrc.org/genomes/${gid}/${gid}.PATRIC.faa
hmmsearch --cpu 4 --tblout ${gid}.tbl -o /dev/null ${profile} ${gid}.PATRIC.faa 
