!#/bin/bash

for genome_accession in "$@"; do 
	wget ftp://ftp.patricbrc.org/genomes/$genome_accession/${genome_accession}.fna
done
