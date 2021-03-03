#!/bin/bash

mash sketch -r ${1}*.gz -o $1 2> /dev/null
mash dist /project/biocomplexity/fungcat/anwarren_locker/mash/patric_all.msh ${1}.msh 2> /dev/null | sort -k3 2> /dev/null | head -n 1 | cut -f 1 | cut -f 6 -d "/" #&& rm ${1}.msh
