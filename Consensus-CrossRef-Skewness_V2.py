#!/usr/bin/env python

"""
This script accepts Kraken/Kaiju output and a file listing all Pearson coefficient 
of skewness, and gives a .tsv file with predictions plus skewness scores
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys
import argparse

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
	k,v = (line[0].strip().split("."))[2], (str(line[1]) + "___" + str(line[2]))
	SkewDict[(k)] = (v)

for prediction in Consensus:
	read, predictedtaxid = prediction[1], prediction[2]
	if predictedtaxid in SkewDict:
		PearsonSkewnessAndHitDist = ((SkewDict[predictedtaxid]).split("___"))
		PearsonSkewness = PearsonSkewnessAndHitDist[0]
		HitDist = PearsonSkewnessAndHitDist[1]
		OutputAll.write("\t".join(prediction) + "\t" + str(PearsonSkewness) + "\t" + str(HitDist) + "\n")
		try:
			HitDist = float(HitDist)
		except:
			continue
		if (type(HitDist) == float):
			if int(HitDist) > 70: #If more than 70% of pseudo-chromosomes have a hit:
				if -0.6 < float(PearsonSkewness) < 0.6: #If the Skewness is between -0.6 and 0.6
					Output.write("\t".join(prediction) + "\t" + str(PearsonSkewness) + "\t" + str(HitDist) + "\n")

OutputAll.close()
Output.close()
