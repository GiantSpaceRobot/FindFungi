#Script to calculate the Pearson's coefficient of skewness for a set of hit distribution data
#Usage: python2.7 Skewness-Calculator.py ContigLengths_Taxid-File Hit-Distribution-File Skewness.Output

import sys
from math import sqrt
import os

Output = open(sys.argv[3], "w")

ContigLenDict = dict()
NumberOfContigs = 0
ContigsWithHits = 0
#Build a dictionary from contig lengths. key = contig, val = length
with open(sys.argv[1]) as ContigLengths:
	for line in ContigLengths:
		NumberOfContigs = NumberOfContigs + 1
		(key, val) = line.strip().split(" ")
		ContigLenDict[(key)] = (val)
#Add hit distribution data to the dictionary
if os.stat(sys.argv[2]).st_size == 0:
	print "The file %s is empty. Cannot calculate Pearson Coefficient of Skewness" % sys.argv[2]
	Taxid = (str(sys.argv[1]).strip().split("-"))[-1]
	Output.write(Taxid + "\tN/A (no BLAST hits)\n")
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
		if len(v.split("\t")) == 1:  #If there is no value for v, make it 0. This indicates that no reads mapped to that contig.
			v = str(v) + "\t0"
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
	if (float(ContigsWithHits)/float(NumberOfContigs)) > 0.7:
		Output.write(Taxid + "\t" + str(PearsonSkewness) + "\n")
	else:
		print "No BLAST hits on at least 70 percent of contigs for %s" % (sys.argv[2])
		Output.write(Taxid + "\t" + str(PearsonSkewness) + " (no BLAST hits on at least 70% of contigs)\n")
Output.close()
