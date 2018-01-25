#!/bin/sh
# Usage: LowestCommonAncestor_V4.sh Results-File.tsv
# Author: Paul Donovan 
# Email: pauldonovandonegal@gmail.com
# 25-Jan-2018

x=$1
y=${x%.*} 
awk '{print $3}' $1 | sort | uniq  > $y-taxids.txt
/home/paul/local/bin/python2.7 /home/paul/scripts/LowestCommonAncestor_V4.py $1 $y-taxids.txt $y-lca.csv
echo "Done"
