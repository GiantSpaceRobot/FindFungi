#!/usr/bin/env python

"""
This script accepts Kraken output and BLAST output, gets the reads 
common to both, and outputs a kraken output file sorted by read
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('KrakenOutput.tsv')
parser.add_argument('BLAST-output-reads')
parser.add_argument('ReducedKrakenOutput.tsv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

#Reads classified as fungal by BLAST are read into a list
BLASTResultsList = list(line.strip().split("\t") for line in open(sys.argv[2]))

#Open the output file for writing
Output = open(sys.argv[3], "w")

#Define dictionaries
ReadDict = dict()
BLASTDict = dict()

#Add BLAST reads as dictionary keys with values of 0
for read in BLASTResultsList:
	BLASTDict[read[0]] = [0]

#Loop over Kraken output, find reads common to Kraken and BLAST, output in Kraken format
for line in open(sys.argv[1], "r"):
	i = line.strip().split("\t")
	read = i[1]
	Line = "\t".join(i)
	if str(read) in BLASTDict:
		Output.write(Line + "\n")
Output.close()
