#!/bin/bash

module load gcc bwa
#This script requires seqtk (https://github.com/lh3/seqtk/releases), and
#RSeQC (infer_experiment.py, can be installed via pip)

#Input is a reference fasta file and the reads from a single or paired end run
#A bed file of same prefix as the fasta reference must exist and can be built from gff file
#This can be done using bedopts gff2bed https://bedops.readthedocs.io/en/latest/content/installation.html#linux

#Fasta formatted reference for alignment with BWA
reference=$1
shift

#The rest of the inputs should either be a single read file or the two pairs of a paired end run
for read_file in "$@"; do
	#Make a subsampled version of the file in the current directory
	file=${read_file##*/}
	seqtk sample -s100 $read_file 1000 > ${file}.subsampled.temp
done

#Align them to the reference
bwa mem -t 4 $reference *.subsampled.temp 2> /dev/null > subsampled_output.bam
rm *.subsampled.temp

infer_experiment.py -i subsampled_output.bam -r ${reference%.f*}.bed -s 20000 2> /dev/null | tail -n 4 | sed 's/+":/+" (fr): /' | sed 's/-":/-" (rf): /'

rm subsampled_output.bam
