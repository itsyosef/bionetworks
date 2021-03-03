#!/bin/bash

#Download the necessary files for the reference sequnce
echo 'Processing e coli k-12 substrain MG1655'

annogesic create --project_path coli

annogesic get_input_files --project_path coli --ftp_path  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/ --ref_gff --ref_gbk

python3 get_lib_line.py coli

annogesic transcript --annotation_files coli/input/references/annotations/NC_000913.3.gff --tex_notex_libs $(python3 get_lib_line.py coli) --project_path coli --replicate_tex 1_1


annogesic operon --annotation_files coli/input/references/annotations/NC_000913.3.gff --transcript_files coli/output/transcripts/gffs/NC_000913.3_transcript.gff  --project_path coli 