#!/usr/bin/env python

"""
Script that takes in .csv file and outputs an .R script to make a wordcloud
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

import sys
import re
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Sorted-LCA.csv')
parser.add_argument('Output.R')
parser.add_argument('/Some/file/location.pdf')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

#The number of reads required to keep a taxid
CutOff = 50  

#Read CSV file into a list
File_Read = File_Read = list(line.strip().split(",") for line in open(sys.argv[1]))

#Open output R file for writing
File_out = open(sys.argv[2], "w")

#Write R script header
File_out.write('library(wordcloud)\n\na=c(')

#Define lists and variables
PDF_location = str(sys.argv[3])
TaxonList = list()
FreqList = list()
BadChars="#$[]{}()|"
Dataset = str((sys.argv[2]).strip(".WordCloud.R")).split("/")[-1]

#Reformat taxonomy names and add them to lists
for line in File_Read:
	if line[0].startswith("#"):
		pass
	else:
		if int(line[2]) > int(CutOff):
			BuildTaxon = line[0]                                        #
			for char in BuildTaxon:										#
				if char in BadChars:									#
					BuildTaxon = BuildTaxon.replace(str(char),"") 		#
			BuildTaxon = BuildTaxon.split(" ")							#
			if len(BuildTaxon) > 1:										# Replace long and ugly species/strain
				if BuildTaxon[1] == "sp.":								# names with shorter cleaner names
					Taxon = (BuildTaxon[0]) + " sp. " + BuildTaxon[2] 	#
				else:													#
					Taxon = (BuildTaxon[0])[0] + ". " + BuildTaxon[1]	#
			else:														#
				Taxon = BuildTaxon[0]									#
			TaxonList.append(Taxon)
			FreqList.append(line[2])

#Write the taxonomy names to the R script
for line in TaxonList:
	taxon = '"' + str(line) + '"'
	if line == TaxonList[-1]: #If this is the last line of the file, do not add another comma
		File_out.write(taxon + ")\n\nb=c(")
	else:
		File_out.write(taxon + ",")

#Write frequencies to the R script and finish the script
for index, line in enumerate(FreqList):
	LastIndex = len(FreqList) - 1 #The index of the last item in the list
	if index == LastIndex:
		File_out.write(str(line + ')\n\npdf("' + PDF_location + '", width=8.5, height=8.5, onefile=FALSE)\nx <- paste("Top group = ", a[1])\nreads <- prettyNum(b[1], big.mark=",", scientific=FALSE)\ny <- paste("Number of reads = ", reads)\npar(bg="black")\npar(fg="white")\nwordcloud(a, b, col=terrain.colors(length(a), alpha=1), rot.per=0, random.order = FALSE)\nlegend("bottomleft", c(x, y), col="white", bty = "n")\ntitle(main = "' + Dataset + '", col.main="white")\ndev.off()'))
	else:
		File_out.write(str(line) + ",")

File_out.close()

