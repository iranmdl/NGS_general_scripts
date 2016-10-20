#!/bin/bash 
file="$1"
cat $file | grep "GO:" | cut -f1,14 | awk -F'\t' '{split($2,a,/\|/); for (x in a){print $1"\t"a[x]}}' | sort -u
