#!/bin/bash

query=$1
query=${query%.*}

grep -P "\t$query" /scratch/jho5ze/bionets/operons/assembly_summary_refseq.txt | cut -f 20 | head -n 1
