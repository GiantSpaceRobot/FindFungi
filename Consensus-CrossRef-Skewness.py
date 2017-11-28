#!/usr/bin/env python

"""
This script accepts Kraken/Kaiju output and a file listing all Pearson coefficient 
of skewness, and gives a .tsv file with predictions plus skewness scores
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Consensus_KrakenOutput.tsv')
parser.add_argument('Skewness-Scores')
parser.add_argument('Output')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

Consensus = list(line.strip().split("\t") for line in open(sys.argv[1]))
Skewness = list(line.strip().split("\t") for line in open(sys.argv[2]))
Output = open(sys.argv[3], "w")
OutputAll = open(str(sys.argv[3]) + "_AllResults.tsv", "w")

SkewDict = dict()

for line in Skewness:
	k,v = line[0], line[1]
	SkewDict[(k)] = (v)

for prediction in Consensus:
	read, predictedtaxid = prediction[1], prediction[2]
	if predictedtaxid in SkewDict:
		PearsonSkewness = SkewDict[predictedtaxid]
		OutputAll.write("\t".join(prediction) + "\t" + str(PearsonSkewness) + "\n")
		if "BLAST" not in str(PearsonSkewness): #Check if the Pearson Skewness variable is a number (by checking if it does not contain the string BLAST)
			if -0.3 < float(PearsonSkewness) < 0.3: #If the Skewness is between -0.3 and 0.3
				Output.write("\t".join(prediction) + "\t" + str(PearsonSkewness) + "\n")

OutputAll.close()
Output.close()
