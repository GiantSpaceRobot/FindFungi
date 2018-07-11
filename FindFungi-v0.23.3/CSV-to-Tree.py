#!/usr/bin/env python

"""
This script accepts .csv pipeline output and gives a .ps file with a basic tree structure
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys
import argparse
from ete3 import NCBITaxa

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Sorted-LCA.csv')
parser.add_argument('Output.gv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

ncbi = NCBITaxa()

#The number of species you want to create the tree with
NumberOfSpecies = 10   

#Read CSV results into list, remove all but the top 10 most abundant taxonomies
ResultsList = list(line.strip().split(",") for line in open(sys.argv[1]))
ResultsList = ResultsList[0:int(NumberOfSpecies) + 1] #Take first n items in list (+1 is to negate the header line)

#Open output file for writing
Output = open(sys.argv[2], "w")

#Write header line in dot format
Output.write('digraph G {\n\tsize="8,5!";\n')

#Define lists, dicts and variables
ResultTaxids = list()
TreeList = list()
BadChars = "()[]{}/|"
TaxidFreqDict = {}
Counter = 0

#Re-open CSV file, create a dictionary with taxid as key and number of reads as value
with open(sys.argv[1]) as f:
	for line in f:
		if not line.startswith("#"):
			tok = line.strip().split(",")
			TaxidFreqDict[tok[1]] = tok[2]

#Build the dot script
for line in ResultsList:
	if line[0].startswith("#"):
		pass
	else:
		ResultTaxid = line[1]
		ResultTaxids.append(ResultTaxid)
		lineage = ncbi.get_lineage(ResultTaxid)
		for index, taxid in enumerate(lineage):
			name = ncbi.get_taxid_translator([str(taxid)])
			name = name[taxid]
			for char in name:
				if char in BadChars:
					name = name.replace(str(char),"_") #Replace ugly strings
			NextIndex = int(index) + 1
			if NextIndex == len(lineage):
				pass
			else:
				NextTaxid = lineage[NextIndex]
				NextName = ncbi.get_taxid_translator([str(NextTaxid)])
				NextName = NextName[NextTaxid]
				for char in NextName:
					if char in BadChars:
						NextName = NextName.replace(str(char),"_") #Replace ugly strings
				NodeToNode = str('\t"' + str(name) + '" -> "' + str(NextName) + '";\n')  
				if any(NodeToNode in s for s in TreeList):
					pass
				else:
					Output.write(NodeToNode)	
					TreeList.append(NodeToNode)
					if str(NextTaxid) in TaxidFreqDict: #If there is information available about number of reads for this taxid, use it
						value = TaxidFreqDict[str(NextTaxid)]
						Freq = format(int(value), ",d")  #Adds commas to numbers to make them more human-readable
						Output.write(str('\t"' + str(NextName) + '" [xlabel="' + str(Freq) + ' reads"];\n'))
			
Output.write("}")
Output.close()
