#!/bin/bash

module load gcc bwa

wget -qN -O ${1}.fna ftp://ftp.patricbrc.org/genomes/${1}/${1}.fna
bwa index ${1}.fna

wget -qN -O ${1}.gff ftp://ftp.patricbrc.org/genomes/${1}/${1}.PATRIC.gff && sed -i 's/^accn|//g' ${1}.gff

awk '{Q=$0;gsub(/CDS/, "gene", Q);gsub(/ID=/, "ID=gene_for-", Q); if ($0 ~ "CDS") print $0 ORS Q; else print $0}' ${1}.gff > ${1}.gff.temp; mv ${1}.gff.temp ${1}.gff

gff2bed < ${1}.gff > ${1}.bed
