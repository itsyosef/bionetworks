#!/bin/bash

rg "hCoV-19/USA" -A 1 $1 > $1.new_usa

rg -v "^--$" $1.new_usa > $1.new_usa.nodash

sed -r "s/>hCoV-19\/USA\//>/g" $1.new_usa.nodash | sed 's/.*|E/>E/g' | tr "|" " " >> $1.new

mv $1.new $1

rm $1.new_usa

rm $1.new_usa.nodash

screed db $1
