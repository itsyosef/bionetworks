#!/bin/bash

for i in $(cat ../notebooks/antismash_accession_ids.txt ); do rg $(rg $i prokaryotes.txt | cut -f 19) ../notebooks/genome_metadata | cut -f 1 >> antismash_genome_from_accession.txt; done
