#!/usr/bin/env python

"""
This script accepts Kraken/Kaiju output and a file listing
all taxids, and gives a .csv file with taxid statistics
"""

__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys
from ete3 import NCBITaxa
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('KrakenOutput.tsv')
parser.add_argument('Taxids.txt')
parser.add_argument('Output.tsv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

ncbi = NCBITaxa()

ResultsTuple = tuple(line.strip().split("\t") for line in open(sys.argv[1]))
AllTaxidSet = tuple(line.strip() for line in open(sys.argv[2]))
Output = open(sys.argv[3], "w")
Output.write("#Taxon name,Taxid,Reads mapping to taxid,Reads mapping to children taxids,Pearson skewness score,Percent of pseudo-chromosomes with read hits\n")

AllTaxidDict = dict.fromkeys(AllTaxidSet, 0)
ResultTaxids = list()
TaxidPearsonDict = dict()

for line in ResultsTuple:
    ResultTaxid = line[2]
    Pearson = line[5] + "___" + line[6]
    ResultTaxids.append(ResultTaxid)
    if ResultTaxid in AllTaxidDict:
        OldValue = int(AllTaxidDict[ResultTaxid])
        AllTaxidDict[ResultTaxid] = OldValue + 1
    if ResultTaxid in TaxidPearsonDict:
        pass
    else:
        TaxidPearsonDict[ResultTaxid] = Pearson
for k,v in AllTaxidDict.iteritems():
    ChildrenCount = 0
    try:
        descendants = ncbi.get_lineage(str(k))
    except ValueError:
        pass
    taxid = list([k])
    name = ncbi.get_taxid_translator(taxid)
    for i in descendants:
        if str(i) in AllTaxidDict:
            ChildrenCount = ChildrenCount + int(AllTaxidDict[str(i)])
    PearsonSkewness = TaxidPearsonDict.get(k)
    Skewness = (PearsonSkewness.strip().split("___"))[0]
    HitDist= (PearsonSkewness.strip().split("___"))[1]
    if PearsonSkewness == "":
        Output.write(str(name.itervalues().next()) + "," + str(k) + "," + str(v) + "," + str(ChildrenCount) + ",N/A (Pearson score only calculated for leaf nodes\n")
    else:
        Output.write(str(name.itervalues().next()) + "," + str(k) + "," + str(v) + "," + str(ChildrenCount) + "," + str(Skewness) + "," + str(HitDist) + "\n")
Output.close()
