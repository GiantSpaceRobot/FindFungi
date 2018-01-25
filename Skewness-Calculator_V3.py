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

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Results.tsv')
parser.add_argument('BLAST-Hit-Distribution.tsv')
parser.add_argument('Output.tsv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

Output = open(sys.argv[3], "w")

ContigLenDict = dict()
NumberOfContigs = 0
ContigsWithHits = 0
#Build a dictionary from contig lengths. key = contig, val = length
with open(sys.argv[1]) as ContigLengths:
	for line in ContigLengths:
		print line
		NumberOfContigs = NumberOfContigs + 1
		(key, val) = line.strip().split("\t")
		print (key, val)
		ContigLenDict[(key)] = (val)
#Add hit distribution data to the dictionary
if os.stat(sys.argv[2]).st_size == 0:
	print "The file %s is empty. Cannot calculate Pearson Coefficient of Skewness" % sys.argv[2]
	Taxid = (str(sys.argv[1]).strip().split("-"))[-1]
	Output.write(Taxid + "\tN/A\t(no BLAST hits)\n")
else:
	with open(sys.argv[2]) as HitDistribution:
		for line in HitDistribution:
			ContigsWithHits = ContigsWithHits + 1
			(v, k) = line.strip().split(" ")
			if k in ContigLenDict:
				NewVal = str(ContigLenDict[k]) + "\t" + str(v)
				ContigLenDict[(k)] = (NewVal)
	ReadsPerContigList = list()
	for k,v in ContigLenDict.iteritems():
		if len(v.split("\t")) == 1:
			v = str(v) + "\t0"   # Uses all contigs, not just those with mapped reads
		(Contig, Reads) = v.strip().split("\t")
		ReadsPerContig = (float(Reads))/(float(Contig))
		ReadsPerContigList.append(ReadsPerContig)
	ReadsPerContigList = sorted(ReadsPerContigList, key=float)
	ListLen = len(ReadsPerContigList)
	if ListLen % 2 == 1:
		Median = sorted(ReadsPerContigList)[ListLen//2]
	else:
		Median = sum(sorted(ReadsPerContigList)[ListLen//2-1:ListLen//2+1])/2.0
	Mean = float(sum(ReadsPerContigList)) / max(len(ReadsPerContigList), 1)
	StandardDev = sqrt(float(reduce(lambda x, y: x + y, map(lambda x: (float(x) - Mean) ** 2, ReadsPerContigList))) / len(ReadsPerContigList))
	PearsonSkewness = (Mean-Median)/StandardDev
	Taxid = (str(sys.argv[1]).strip().split("-"))[-1]
	PercentChroms = (float(ContigsWithHits)/float(NumberOfContigs)*100)
	Output.write(Taxid + "\t" + str(PearsonSkewness) + "\t" + str(PercentChroms) + "\n")

Output.close()
