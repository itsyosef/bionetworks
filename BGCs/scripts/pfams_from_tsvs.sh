#!/bin/bash

# $1=first_pass_split_0_train_with_pfams.tsv
# $2=first_pass_split_0_test_with_pfams.tsv

rg_str="PF"

for gid in $(rg $rg_str $1 | cut -f 1 | sort -u | cut -f 2 -d '"'); do echo $(rg $gid $1 | rg $rg_str | cut -f 18 | tr -d '"' | tr "\n" " ") > ${gid}_pfams.txt; done

for gid in $(rg $rg_str $2 | cut -f 1 | sort -u | cut -f 2 -d '"'); do echo $(rg $gid $2 | rg $rg_str | cut -f 18 | tr -d '"' | tr "\n" " ") > ${gid}_pfams.txt; done
