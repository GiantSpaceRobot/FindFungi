#!/bin/sh
# Usage: FindFungi_GoButton_v0.23.sh
# Author: Paul Donovan 
# Email: pauldonovandonegal@gmail.com
# 25-Oct-2017

echo "Started at $(date)"

#Set working directory and enter it
DIR="/home/paul/Project-work/Metagenomics-Datasets/EBI-Metagenomics"
cd $DIR

#Shuffle queue list, open a file, generate metadata file, download dataset
for d in `find $DIR/queue/ -maxdepth 1 -type f | shuf`; do	
	File=$(basename $d)	
	URL=$(awk -F "\t" '{print $1}' $d)
	URL_0="${URL}${File}.fastq.gz"
	URL_1="${URL}${File}_1.fastq.gz"
	URL_2="${URL}${File}_2.fastq.gz"
	Data=$(awk -F "\t" '{print $2}' $d)
	Description=$(awk -F "\t" '{print $3}' $d)
	Type=$(awk -F "\t" '{print $4}' $d)
	Sequencer=$(awk -F "\t" '{print $5}' $d)
	echo "Working on $File"
	echo "Dataset name: "$Data >> $DIR/Metadata.$File.txt
	echo "Description: "$Description >> $DIR/Metadata.$File.txt
	echo "Data type: "$Type >> $DIR/Metadata.$File.txt
	echo -e "Sequencing platform: "$Sequencer"\n" >> $DIR/Metadata.$File.txt
	if [[ "$File" =~ ^mgm.* ]]; then  #If file is from MG-RAST
		mkdir $File
		mg-download.py --metagenome $URL --dir $File #Use MG-RAST bulk download script
		Downloads=$(find $File -name "*.050.upload.fastq*" -type f)
		echo -e "Download files: "$Downloads
		if [[ "$Downloads" == 0 ]]; then    #If there is no fastq file, go to next iteration
			echo "Error: Skipping "$File", no FASTQ file"
			continue
		elif [[ "$Downloads" -gt 1 ]]; then
			echo "Error: Skipping "$File", more than one FASTQ file"
			continue
		else
			for m in $Downloads; do
				FASTQ=$(basename $m)
				mv $m $DIR/download-data/$FASTQ.$File
			done
		rm -rf $File	
		fi
	else  #If the file is not from MG-RAST, it is an EBI Metagenomics file
		wget -t 2 -c -o $DIR/$File.wget.log0 --passive-ftp $URL_0 & #http://www.linuxquestions.org/questions/solaris-opensolaris-20/wget-seems-to-have-problem-with-ftp-servers-520840/
		pid0=$!
		wget -t 2 -c -o $DIR/$File.wget.log1 --passive-ftp $URL_1 & 
		pid1=$!
		wget -t 2 -c -o $DIR/$File.wget.log2 --passive-ftp $URL_2 &
		pid2=$!
		wait $pid0 $pid1 $pid2 #Attempt to download files, wait for completion
		Downloads=$(find . -maxdepth 1 -name "*fastq*" -type f)
		for m in $Downloads; do
			FASTQ=$(basename $m)
			mv $m $DIR/download-data/
		done
		ssh in01 bsub -K -q C gunzip $DIR/download-data/*$File* #Gunzip the FASTQ file(s)
	fi	
	cat $DIR/download-data/*$File* > $DIR/download-data/$File.Final.fastq #Concatenate FASTQ files
	while [ $(ssh in01 qstat -l | grep Request | grep "bin/kraken" | wc -l) -gt 0 ] 	########
	do																					#
		echo "Waiting for Kraken to finish..."											# Wait for all Kraken 
		sleep 120																		# processes to finish
	done																				#
	echo "Kraken is not running. While loop to see if there are more than 64 jobs"
	while [ $(ssh in01 qstat -l | grep Request | wc -l) -gt 64 ] 						########
	do																					#  
		echo "There are over 64 jobs, waiting for a few to finish before proceeding..."	# If there are many bsub
		sleep 120																		# jobs running, wait
	done																				#
	echo "There are less than 64 jobs. Starting analysis"
	echo "Beginning analysis of "$File >> $DIR/StartingAnalysis.txt
	sendmail email@address.com < $DIR/StartingAnalysis.txt
	rm $DIR/StartingAnalysis.txt
	ssh in01 /home/paul/scripts/./FindFungi-0.23.sh $DIR/download-data/$File.Final.fastq $File & #Begin FindFungi pipeline
	mv $d $DIR/finished/
	while [ ! -d /home/paul/Project-work/Analysed_Datasets/$File/FASTA ]	#
	do																		# Wait for the FASTA directory to be 
		sleep 10															# created before removing the download data
	done																	# 
	rm -f $DIR/download-data/*$File* $DIR/download-data/*prelim*
done

echo "Finished at $(date)"

