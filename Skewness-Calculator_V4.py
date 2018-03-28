#!/usr/bin/env python

"""
Script to calculate the Pearson's coefficient 
of skewness for a set of hit distribution data
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys
from math import sqrt
import os
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Hit-Distribution-Results.tsv')
parser.add_argument('SkewnessScore.tsv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

Output = open(sys.argv[2], "w")

NumberOfContigs = 20
ContigsWithHits = 0

ContigHitDict = dict()

#Add hit distribution data to the dictionary
if os.stat(sys.argv[1]).st_size == 0:
	print "The file %s is empty. Cannot calculate Pearson Coefficient of Skewness" % sys.argv[2]
	Taxid = (str(sys.argv[1]).strip().split("-"))[-1]
	Output.write(Taxid + "\tN/A\t(no BLAST hits)\n")
else:
	with open(sys.argv[1]) as HitDistribution:
		for line in HitDistribution:
			line1 = line.strip().split(" ")
			contig, contigHits = line1[1], float(line1[0])
			ContigHitDict[(contig)] = (contigHits)
			ContigsWithHits = ContigsWithHits + 1
	ReadsPerContigList = list()
	contigsList = list()
	for k,v in ContigHitDict.iteritems():
		contig = k.split('Fragment-')
		contigsList.append(int(contig[-1]))
	for i in range(1,21):
		if int(i) in contigsList:
			pass
		else:
			EmptyContig = contig[0] + "Fragment-" + str(i)
			EmptyContigValue = float(0)
			ContigHitDict[EmptyContig] = EmptyContigValue		
	ReadsPerContigList = ContigHitDict.values()
	[float(i) for i in ReadsPerContigList]
	ReadsPerContigList = sorted(ReadsPerContigList, key=float)
	ListLen = len(ReadsPerContigList)
	if ListLen % 2 == 1:
		Median = sorted(ReadsPerContigList)[ListLen//2]
	else:
		Median = sum(sorted(ReadsPerContigList)[ListLen//2-1:ListLen//2+1])/2.0
	Mean = float(sum(ReadsPerContigList)) / max(len(ReadsPerContigList), 1)
	StandardDev = sqrt(float(reduce(lambda x, y: x + y, map(lambda x: (float(x) - Mean) ** 2, ReadsPerContigList))) / len(ReadsPerContigList))
	PearsonSkewness = 3*(Mean-Median)/StandardDev
	Taxid = (str(sys.argv[1]).strip().split("-"))[-1]
	PercentChroms = (float(ContigsWithHits)/float(NumberOfContigs)*100)
	Output.write(Taxid + "\t" + str(PearsonSkewness) + "\t" + str(PercentChroms) + "\n")
Output.close()
