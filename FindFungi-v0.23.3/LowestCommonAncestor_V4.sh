#!/bin/bash
# Usage: LowestCommonAncestor_V4.sh Results-File.tsv
# Author: Paul Donovan 
# Email: pauldonovandonegal@gmail.com
# 11-Jul-2018

##### USER INPUT REQUIRED:
### Edit this ScriptPath to point to the directory containing the downloaded scripts
ScriptPath=/home/user/scripts  #Location of downloaded python and shell scripts

x=$1
y=${x%.*} 
awk '{print $3}' $1 | sort | uniq  > $y-taxids.txt
python2.7 $ScriptPath/LowestCommonAncestor_V4.py $1 $y-taxids.txt $y-lca.csv
echo "Done"
