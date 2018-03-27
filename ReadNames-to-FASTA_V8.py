#!/usr/bin/env python

"""
Script that takes in read names and a FASTA file, and 
outputs the FASTA sequences corresponding to the reads
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

from Bio import SeqIO
import sys
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('ReadNames.txt')
parser.add_argument('FASTA-input.fsa')
parser.add_argument('FASTA-output.fsa')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

Output = open(sys.argv[3], 'w')

newDict = dict()
with open(sys.argv[1], 'r') as f:
	for line in f:
		if line.startswith('#'):
			pass
		else:
			item = line.strip().split("\t")
			newDict[item[1]] = (item[2])

for fasta in SeqIO.parse(open(sys.argv[2]),'fasta'):
	name, sequence = fasta.id, str(fasta.seq)
	if name in newDict:
		Output.write(">" + str(name) + "\tTaxon-prediction: " + str(newDict.get(name)) + "\n" + str(sequence.split(":")[0]) + "\n")
Output.close()
