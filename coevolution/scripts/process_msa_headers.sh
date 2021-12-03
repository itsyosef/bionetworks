#!/bin/bash

sed -r "s/>hCoV-19\/USA\//>/g" $1 | sed 's/.*|E/>E/g' | tr "|" " " >> $1.new

mv $1.new $1

screed db $1
