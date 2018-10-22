
# FindFungi-v0.23.3

A pipeline for the identification of fungi in public metagenomics datasets.
The FindFungi pipeline uses the metagenomics read-classifier Kraken with 32 custom fungal 
databases to generate 32 taxon predictions for a single read. These 32 predictions are 
combined to generate a consensus prediction. All reads are then BLASTed against their
predicted genomes to generate read distribution skewness scores to select for the 
most likely true positives.

FindFungi-v0.23.3 corrects an error in the LowestCommonAncestor_V4.sh script where the 
Python path and script path were incorrect.

FindFungi v0.23.2 corrects an error where the KrakenReduction.py script requested a file 
that was generated in FindFungi v0.22. This step has been removed.

FindFungi v0.23.1 corrects an error in v0.23, which did not correctly calculate Pearson’s 
skewness scores. A multiplication factor of 3 was omitted in v0.23 
(3(mean-median)/standard deviation). Skewness scores calculated by v0.23.1 are therefore 
3-times higher than with v0.23 (i.e. cut-offs should be expressed as -0.6 to +0.6, and 
not -0.2 to +0.2 as in v0.23). This makes no difference to the interpretation of any of 
the results or the data analysis. 

FindFungi-0.23 was built on an IBM platform load-sharing facility with 32 worker nodes.
Similar architecture is required to set up the pipeline due to the memory requirements of 
the Kraken databases.

## Quickstart

Download the pipeline, databases, associated scripts, prerequisites and other tools. 
Run the pipeline:

```
./FindFungi-0.23.3.sh /path/to/FASTQ-file.fastq Dataset-name
```

## Getting Started

These instructions will hopefully allow you to get a copy of FindFungi up and running 
on your own compute-cluster/server for development or your own analyses. If using a 
non-IBM LSF compute cluster, change the 'bsub' commands to reflect your architecture.
If using a single server, remove the 'bsub' commands.

### Prerequisites

FindFungi v0.23 was built using the following:

* gcc version 4.4.4 20100726 (Red Hat 4.4.4-13)
* coreutils 8.27
* python 2.7.13 (modules: sys, os, ete3, Bio, math, argparse, itertools, collections, 
re)
* skewer 0.2.2
* kraken 0.10.5-beta
* ncbi blast 2.2.30
* Rscript 3.3.3 (packages: wordcloud)
* graphviz 2.40.1

### Installing

* Download all of the scripts from GitHub/GiantSpaceRobot and move to a directory 
(/your/directory/scripts). You may need to give these scripts more permissions (e.g. chmod 755 *).
* In the FindFungi-v0.23.3 script, change the absolute paths of skewer, kraken, blast, the shell 
and python scripts to reflect your environment, or add these tools and scripts to you $PATH. 
You will also need to edit the LowestCommonAncestor.sh script to include the path to the downloaded scripts.
* NOTE: It may be necessary for you to include the absolute paths for all of the scripts
and tools within the FindFungi-0.23.sh master script, depending on the cluster node
preferences (e.g. executing 'python' actually calls the node's version of python, not yours).
* Download the Kraken and BLAST databases from this website (http://bioinformatics.czc.hokudai.ac.jp/findfungi/).
* Uncompress these files and put them somewhere sensible:

```
tar -xvfz Kraken_*.tar.gz
mv Kraken_* Kraken_DB_Directory/
```

### Testing the pipeline

Download the dataset ERR675624 from the European Nucleotide Archive database. This 
dataset contains fungal reads.

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR675/ERR675624/ERR675624_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR675/ERR675624/ERR675624_2.fastq.gz
```

Gunzip these files and concatenate them, as we have no need to read pair information.

```
gunzip ERR675624_*.fastq.gz
cat ERR675624_*.fastq > ERR675624_both-pairs.fastq
```

Execute the FindFungi pipeline on this FASTQ file. We use nohup here to allow the 
pipeline to run in the background.

```
nohup ./FindFungi-0.23.sh /path/to/ERR675624_both-pairs.fastq ERR675624
```

The first command line argument (/path/to/ERR675624_both-pairs.fastq) points to your FASTQ
file. The second (ERR675624) is the name FindFungi will use for this dataset. This name 
should be informative.

The .csv results should show the following:

```
#Taxon name,Taxid,Reads mapping to taxid,Reads mapping to children taxids,Pearson skewness score,Percent of pseudo-chromosomes with read hits
Candida sp. LDI48194,1759314,671,0,0.524623062587,100.0
Malassezia restricta,76775,378,0,0.496034792692,100.0
Candida tropicalis MYA-3404,294747,265,0,-0.265788716977,100.0
```

## Contributors

* Paul Donovan, PhD (email: pauldonovandonegal@gmail.com)
* Gabriel Gonzalez, PhD (e-mail: gagonzalez@czc.hokudai.ac.jp)

## License

This project is licensed under the MIT License.
