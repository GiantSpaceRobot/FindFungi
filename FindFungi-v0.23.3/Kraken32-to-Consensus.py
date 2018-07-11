#!/usr/bin/env python

"""
Script that takes the concatenated sorted output of the Kraken32 
databases and finds a consensus prediction for each read
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

from itertools import islice
import sys
from collections import Counter
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Sorted.Kraken.tsv')
parser.add_argument('Consensus.Kraken.tsv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

#Define lists and variables
n = 32 #The number of lines to loop over at once
GuessList = list()
YesCount = 0
NoCount = 0

#Open output file for writing
Output = open(sys.argv[2], "w")

#Find the consensus prediction using 32 predictions for the same read 
with open(sys.argv[1]) as f:  #http://stackoverflow.com/questions/6335839/python-how-to-read-n-number-of-lines-at-a-time
	while True:
		next_n_lines = list(islice(f, n))  #Loop over 32 lines at a time because the file is often too large to read in
		if not next_n_lines: #If the list is empty, exit loop because we are at the end of the file
			break
		for i in next_n_lines: #Looping over 32 predictions for the same read
			Guess = (i.strip().split("\t"))[2] 
			if str(Guess) == str(0): #If there was no prediction for the read
				pass
			else: #There is a prediction 
				GuessList.append(Guess) #Add all taxid predictions to this list 
		if len(GuessList) == 0: #If there as no prediction at all made for the current read, skip this list of 32
			NoCount = NoCount + 1
			continue       #Move onto the next 32 lines
		else:
			YesCount = YesCount + 1
			count = Counter(GuessList)               #                                       #
			freq_list = count.values()               #                                       #
			max_cnt = max(freq_list)                 # Find the most common item in the list #
			total = freq_list.count(max_cnt)         #                                       #
			most_common = count.most_common(total)   #                                       #
			if len(most_common) == 1:   #If there is one most common taxid
				read = (next_n_lines[0].strip().split("\t")[1])
				Output.write("C\t" + read + "\t" + (most_common[0])[0] + "\t100\t100\n")
			else:  #If there is more than one most common prediction
				PredictionList = list()
				for i in next_n_lines:  #Loop over the 32 kraken predictions again
					item = i.strip().split("\t")
					if str(item[2]) == str(0): #If kraken failed to predict a taxid for the read (aka 0), skip this 
						pass
					else:
						Prediction = item[4].split(" ") #Gather all predictions from all kraken runs
						for u in Prediction:
							if "A" not in str(u): #Don't use the 'ambiguous' flag in our taxonomy predictions
								PredictionList.append(u) #Add all predictions to a list
				NewDict = dict()
				for kmer in PredictionList:												#
					Taxid = (kmer.split(":"))[0]										#
					KmerNumber = (kmer.split(":")[1])									#
					if str(Taxid) == str(0):											#
						pass															# Loop over all predictions from all 32 kraken runs
					else:																# and create a dictionary where taxid is the key 
						if str(Taxid) in NewDict:										# and the number of kmers supporting the taxid is 
							NewDict[Taxid] = (int(KmerNumber) + int(NewDict[Taxid]))	# the value
						else:															#
							NewDict[Taxid] = int(KmerNumber)							#
				if bool(NewDict) == False: #If the dictionary is empty, don't do anything. This should never be the case
					pass			
				else:
					TaxidPrediction = str(max(NewDict, key=NewDict.get)) #Get the key with the highest value. Draws may occur but a key will be selected
					TypicalRead = next_n_lines[0].strip().split("\t")
					Output.write("C\t" + TypicalRead[1] + "\t" + TaxidPrediction + "\t100\t100\n") 
				NewDict.clear()
			del GuessList[:]

print "Number of reads with consensus prediction: %s" % YesCount
print "Number of reads without a prediction: %s" % NoCount

Output.close()
