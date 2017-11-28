#!/bin/sh
# Usage: LowestCommonAncestor_V3.sh Results-File.tsv
# Author: Paul Donovan 
# Email: pauldonovandonegal@gmail.com
# 28-Nov-2017

x=$1
y=${x%.*} 
awk '{print $3}' $1 | sort | uniq  > $y-taxids.txt
/home/paul/local/bin/python2.7 /home/paul/scripts/LowestCommonAncestor_V3.py $1 $y-taxids.txt $y-lca.csv
echo "Done"
