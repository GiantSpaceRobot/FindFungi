# FindFungi

A pipeline for the identification of fungi in public metagenomics datasets.
The FindFungi pipeline uses the metagenomics read-classifier Kraken with 32 custom fungal 
databases to generate 32 taxon predictions for a single read. These 32 predictions are 
combined to generate a consensus prediction. All reads are then BLASTed against their
predicted genomes to generate read distribution skewness scores to select for the 
most likely true positives.

## Quickstart

Download the pipeline, associated scripts, prerequisites and other tools. 
Run the pipeline:

```
./FindFungi-0.23.sh /path/to/FASTQ-file.fastq Dataset-name
```

## Getting Started

These instructions will hopefully allow you to get a copy of FindFungi up and running 
on your own compute-cluster for development or your own analyses.

### Prerequisites

FindFungi v0.23 was built using the following:

*gcc version 4.4.4 20100726 (Red Hat 4.4.4-13)
*python 2.7.13 (modules: sys, os, NCBITaxa, Bio, math, argparse, itertools, collections, 
re)
*skewer 0.2.2
*kraken 0.10.5-beta
*ncbi blast 2.2.30
*Rscript 3.2.2

### Installing

* Download all of the scripts from GitHub/GiantSpaceRobot and move to a directory 
(/your/directory/scripts).
* Change the absolute paths of kraken, blast, the shell and python scripts to reflect
your environment, or add these tools and scripts to you $PATH.

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
Candida sp. LDI48194,1759314,671,0,0.17486834948,100.0
Malassezia restricta,76775,378,0,0.16534366649,100.0
Candida tropicalis MYA-3404,294747,265,0,-0.0886001807602,100.0
```

## License

This project is licensed under the MIT License.
