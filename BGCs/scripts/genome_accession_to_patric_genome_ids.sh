#!/bin/bash

#for i in $(cat ../notebooks/antismash_accession_ids.txt ); do rg $(rg $i prokaryotes.txt | cut -f 19) ../notebooks/genome_metadata | cut -f 1 >> antismash_genome_from_accession.txt; done

misc_data=/sfs/lustre/bahamut/scratch/jho5ze/bionets/BGCs/misc_data

for accession in "$@"; do
	rg $(rg $accession ${misc_data}/prokaryotes.txt | cut -f 19) ${misc_data}/genome_metadata  | rg -v "\tPlasmid" | cut -f 1;
done
