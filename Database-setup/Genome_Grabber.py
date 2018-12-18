##### Script for looping over NCBI fungi directory and finding all representative genomes. 
##### These are copied to a different location, where the FASTA identifiers and filenames are modified to include Taxid and Species name

import subprocess
import os
import sys
from Bio import SeqIO
from Bio import Seq

OriginalWD = os.getcwd()  #The working directory is stored in a variable

Output = open("Output.txt", "w")
Output2 = open("ErrorReporting.txt", "w")

Count = 0

with open("assembly_summary.txt", "r") as Assembly: #Read in the NCBI file including all information about genomes
	for line in Assembly:
		if line.startswith("#"):
			pass
		else:
			if ("representative" in line) or ("reference" in line):  #We only want the genomes belonging to representative species
				Count = Count + 1
				Output.write(line) #Output representative species information to a file
				lineSplit = line.strip().split("\t")
				Output2.write("Number of representative species: " + str(Count) + "\tSpecies: " + lineSplit[7])
				Species = lineSplit[7]
				if "[" in Species:
					Species = Species.replace("[","")
					Species = Species.replace("]","")
				SpeciesNoSpaces = Species.replace(" ", "_")
				ShortSpeciesNoSpaces = (SpeciesNoSpaces.split("_"))[0] + "_" + (SpeciesNoSpaces.split("_"))[1] 
				os.chdir(OriginalWD)  #Make sure we begin in the right directory
				try:    #This section of try/excepts deal with finding the directory belonging to the species name and moving into it
					Directory = str(os.getcwd()) + "/fungi/" + str(SpeciesNoSpaces)
					os.chdir(Directory)
				except:
					try:
						Directory = str(os.getcwd()) + "/fungi/" + str(ShortSpeciesNoSpaces)
						os.chdir(Directory)
					except:
						try:
							Directory = str(os.getcwd()) + "/fungi/_" + str(SpeciesNoSpaces)
							os.chdir(Directory)
						except:
							try:
								Directory = str(os.getcwd()) + "/fungi/_" + str(ShortSpeciesNoSpaces)
								os.chdir(Directory)
							except:
								Output2.write("No directory for %s" % SpeciesNoSpaces)
								continue  #Continue with next iteration of loop
				try: #Enter the 'representative' directory if present
					os.chdir("representative")
				except:
					try:
						os.chdir("reference")
					except:
						Output2.write("No representative or reference directory in %s" % str(os.getcwd()))
						continue
				FindNames = subprocess.check_output(['find', '.', '-name', '*_genomic.fna.gz*']) #Find the path for the representative genome(s)
				FindNames = FindNames.split("\n")
				for item in FindNames:
					if ".fna" not in item:  #Skip anything that isn't a fASTA file
						pass
					else:
						if "from_genomic" in item: #If the file is a cds or rna file, skip it
							pass
						else:
							Taxid = lineSplit[5]
							GenomePath = item[1:] #Path to genome from current directory
							FullPath = str(os.getcwd()) + str(GenomePath) #Full path to genome	
							TheGenome = GenomePath.split("/")[-1] #Genome file (not named after species)
							TheGenomeUnzipped = TheGenome[:-3] #Remove .gz from genome name
							subprocess.check_output(['gunzip', '%s' % FullPath]) #Unzip the genome file so we can edit it
							GunzipFullPath = FullPath[:-3] #Remove the .gz from the variable FullPath 
							FASTAGenome = "__" + TheGenome.replace("_genomic.fna.gz","") + "__|kraken:taxid|" + str(Taxid)  #New string to replace the ID in FASTA file
							try: #Attempt to open the unzipped genome, edit the identifier and write it out in the new location
								with open(GunzipFullPath, "r") as fin:
									with open("/home/paul/Project-work/Databases/Fungal-Genomes/NCBI-Fungi/RepresentativeSpecies/" + SpeciesNoSpaces + "__Taxid-" + str(Taxid) + "__" + TheGenomeUnzipped, "w") as fout:
										for line in fin:
											if ">" in line:
												ID = line.split(" ")
												NewID = ID[0] + FASTAGenome + " " + " ".join(ID[1:])
												fout.write(NewID)
											else:
												fout.write(line)
							except:
								Output2.write("There is a problem opening the input %s \n Or there is a problem writing the output %s" % (GunzipFullPath, str(SpeciesNoSpaces + "__Taxid-" + str(Taxid) + "__" + TheGenomeUnzipped)))
							subprocess.check_output(['gzip', GunzipFullPath]) #Re-zip the original file
Output.close()
Output2.close()
