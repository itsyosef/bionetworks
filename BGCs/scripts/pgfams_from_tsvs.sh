#!/bin/bash

# $1=first_pass_split_0_train.tsv
# $2=first_pass_split_0_test.tsv

for gid in $(rg PGF $1 | cut -f 1 | sort -u | cut -f 2 -d '"'); do echo $(rg $gid $1 | rg PGF | cut -f 18 | cut -f 2 -d '"' | tr "\n" " ") > ${gid}_pgfams.txt; done

for gid in $(rg PGF $2 | cut -f 1 | sort -u | cut -f 2 -d '"'); do echo $(rg $gid $2 | rg PGF | cut -f 18 | cut -f 2 -d '"' | tr "\n" " ") > ${gid}_pgfams.txt; done