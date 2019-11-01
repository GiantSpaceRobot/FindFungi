#!/bin/bash
# Usage: ./FindFungi-0.23.sh /path/to/FASTA-file.fastq File-name
# Author: Paul Donovan 
# Email: pauldonovandonegal@gmail.com
# 11-Jul-2018

echo "Started at $(date)"
StartTime="Pipeline initiated at $(date)" 

### Precautionary measures before running analysis
if [ $# -eq 0 ]; then
	echo "No arguments provided"
	exit 1
fi

###Declare paths for FindFungi analysis, scripts, tools etc.
x=$1
z=$2


##### USER INPUT REQUIRED:
ScriptPath=/home/user/scripts  #Location of downloaded python and shell scripts
PreDir=/home/user/Some_Directory/$z  #Location you want FindFungi to build the analysis
KrakenDir=/home/user/Location_of_Kraken_Databases  #Location of the 32 downloaded Kraken databases
FungTaxDir=/home/user/Location_of_Fungal-Taxids.txt  #Location of the Fungal taxids and PipelineSummary files from GitHub
BLAST_DB_Dir=/home/user/Location_of_BLAST_Databases  #Location of the 949 downloaded BLAST databases
#####


Dir=$PreDir/FindFungi 

if [ ! -d $PreDir ]; then
	mkdir $PreDir
	echo "Number of reads in the raw input: " >> $PreDir/Run_Statistics.txt
	LinesInReadsIn=$(wc -l $x | awk '{print $1}') 
	ReadsIn=$((LinesInReadsIn/4))
	echo $ReadsIn >> $PreDir/Run_Statistics.txt
	mkdir $PreDir/ReadTrimming
	bsub -K -q C skewer -l 30 -q 15 -t 30 -o $PreDir/ReadTrimming/$z $x &
	wait
	echo "Number of reads after trimming: " >> $PreDir/Run_Statistics.txt
	LinesInReadsLeft=$(wc -l $PreDir/ReadTrimming/$z-trimmed.fastq | awk '{print $1}') 
	ReadsLeft=$((LinesInReadsLeft/4))
	echo $ReadsLeft >> $PreDir/Run_Statistics.txt
	mkdir $PreDir/FASTA
	sed -n '1~4s/^@/>/p;2~4p' $PreDir/ReadTrimming/$z-trimmed.fastq > $PreDir/FASTA/$z.fna #Convert FASTQ to FASTA
	LineCt=$(wc -l $PreDir/FASTA/$z.fna | awk '{print $1}')
	SplitN=$((LineCt/32 + 1))
	SplitI=$(printf "%.0f" $SplitN)
	split -l $SplitI $PreDir/FASTA/$z.fna $PreDir/FASTA/Split.
	for d in $PreDir/FASTA/*Split.*; do 
		bsub -K -q C sed -i 's/\ /_/g' $d & #Replace whitespace with underscore
	done
	wait
	cat $PreDir/FASTA/*Split.* > $PreDir/FASTA/$z.final.fna  
	mkdir $Dir
	mkdir $Dir/Processing
	mkdir $Dir/Processing/SplitFiles_Kraken
	mkdir $Dir/Results
	mkdir $Dir/Results/BLAST_Processing
	mkdir $Dir/bsub_reports

### Release the Kraken 
for i in $(seq 1 32); do
	bsub -K -q C kraken --db $KrakenDir/Kraken_$i --threads 30 --fasta-input $PreDir/FASTA/$z.final.fna --output $Dir/Processing/SplitFiles_Kraken/$z.$i &
done	
wait
for d in $Dir/Processing/SplitFiles_Kraken/*; do
    File=$(basename $d)
    grep ^C $d > $Dir/Processing/$File.Classified.tsv
    wait
done
wait 
cat $Dir/Processing/*Classified* | awk '{print $2}' | sort | uniq > $Dir/Processing/AllClassified_$z #Get all classified reads (no duplicates)

### Gather reads predicted as fungal by Kraken and reformat for later use
sed 's/ /_/g' $PreDir/FASTA/$z.fna | awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' | awk -v reads="$Dir/Processing/AllClassified_$z" -F "\t" 'BEGIN{while((getline k < reads)>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > $Dir/Processing/Reads-From-Kraken-Output.$z.fsa
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' $Dir/Processing/Reads-From-Kraken-Output.$z.fsa > $Dir/Processing/Reads-From-Kraken-Output.$z.Reformatted.fsa &

### Determine the number of reads classified by Kraken and divide by 32
LineCount=$(wc -l $Dir/Processing/AllClassified_$z | awk '{print $1}')
SplitNum=$((LineCount/32 + 1))
SplitInt=$(printf "%.0f" $SplitNum)

for d in $Dir/Processing/SplitFiles_Kraken/*; do #Sort individual Kraken output files
	File=$(basename $d)
	#bsub -K -q C sort_parallel --parallel 16 -o $Dir/Processing/SplitFiles_Kraken/sorted_$File -k2,2 $d & 
	sort -k2,2 $d > $Dir/Processing/SplitFiles_Kraken/sorted_$File &

done
wait
#bsub -K -q C sort_parallel --parallel 16 -o $Dir/Processing/sorted.$z.All-Kraken-Results.tsv -m -k2,2 $Dir/Processing/SplitFiles_Kraken/*sorted* & #Merge and sort all Kraken output files
cat $Dir/Processing/SplitFiles_Kraken/*sorted* | sort -k2,2 > $Dir/Processing/sorted.$z.All-Kraken-Results.tsv &
wait
bsub -K -q C python2.7 $ScriptPath/Kraken32-to-Consensus.py $Dir/Processing/sorted.$z.All-Kraken-Results.tsv $Dir/Processing/Consensus.sorted.$z.All-Kraken-Results.tsv &
wait

### Count number of predictions for each taxonomic unit and sort
awk '{print $3}' $Dir/Processing/Consensus.sorted.$z.All-Kraken-Results.tsv | sort -n | uniq -c | sort -k1,1nr > $Dir/Processing/Consensus.sorted.$z.All-Kraken-Results.TopTaxonomies.tsv

### Using Kraken taxid predictions and their respective read-counts,
### BLAST reads against their predicted genomes (cut-off is  10 reads predicted per species)
### Ignore all higher taxa e.g. genus, family, class etc.
while read p; do
    ReadNumber=$(echo $p | awk -F ' ' '{print $1}')
    Taxid=$(echo $p | awk -F ' ' '{print $2}')
    if [ $ReadNumber -ge 10 ]
    then
        if grep -Fxq $Taxid $FungTaxDir/Fungal-Taxids.txt
            then #If the predicted taxid is one of the 949 fungal species
                echo "BLASTing $Taxid"
                #Gather read names for each taxid
                awk -v taxid="$Taxid" '$3 == taxid {print $2}' $Dir/Processing/Consensus.sorted.$z.All-Kraken-Results.tsv > $Dir/Processing/ReadNames.$Taxid.txt &
            else
                echo "Not BLASTing $Taxid, wrong taxonomy level"
        fi
    else
        echo "Not BLASTing $Taxid, too few reads"
    fi
done < $Dir/Processing/Consensus.sorted.$z.All-Kraken-Results.TopTaxonomies.tsv
wait

### Gather FASTA reads for each taxid predicted
for d in $Dir/Processing/ReadNames.*; do
	File=$(basename $d)
	Taxid=$(echo $File | awk -F '.' '{print $2}')
	bsub -K -q C -e $PreDir/FindFungi/bsub_reports/ReadNames-to-FASTA.$Taxid.stderr -o $Dir/Processing/ReadNames_bsub.$Taxid.fsa awk -v reads="$Dir/Processing/ReadNames.$Taxid.txt" -F "\t" 'BEGIN{while((getline k < reads)>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' $Dir/Processing/Reads-From-Kraken-Output.$z.Reformatted.fsa &
done
wait

### BLAST against the genome of the predicted species
for d in $Dir/Processing/ReadNames_bsub.*.fsa; do
	File=$(basename $d)
	Taxid=$(echo $File | awk -F '.' '{print $2}')
	tail -n +31 $d | head -n -6 > $Dir/Processing/ReadNames.$Taxid.fsa
	bsub -K -q C blastn -task megablast -query $Dir/Processing/ReadNames.$Taxid.fsa -db $BLAST_DB_Dir/Taxid-$Taxid -out $Dir/Results/BLAST_Processing/BLAST.$Taxid -evalue 1E-20 -num_threads 30 -outfmt 6 &
done
wait

### Get the best hit for every read and determine how many chromosomes these reads hit, and how many times
### Calculate the Pearson coefficient of skewness for each species, and gather together all species skewness scores
for d in $Dir/Results/BLAST_Processing/BLAST*; do
	File=$(basename $d)
	Taxid="${File#BLAST.}"
	awk '! a[$1]++' $d > $Dir/Results/BLAST_Processing/Top-Hits.$File
	awk '{print $2}' $Dir/Results/BLAST_Processing/Top-Hits.$File | sort | uniq -c > $Dir/Results/BLAST_Processing/Hit-Distribution.$File
	bsub -K -q C python2.7 $ScriptPath/Skewness-Calculator_V4.py $Dir/Results/BLAST_Processing/Hit-Distribution.$File $Dir/Results/BLAST_Processing/Skewness.$File &
done
wait
cat $Dir/Results/BLAST_Processing/Skewness* > $Dir/Results/BLAST_Processing/All-Skewness-Scores &

### Combine Kraken results with Skewness scores
bsub -K -q C python2.7 $ScriptPath/Consensus-CrossRef-Skewness_V2.py $Dir/Processing/Consensus.sorted.$z.All-Kraken-Results.tsv $Dir/Results/BLAST_Processing/All-Skewness-Scores $Dir/Results/Final_Results_$z.tsv & 
wait

### Gather all taxonomical predictions and reformat to parsable format
$ScriptPath/LowestCommonAncestor_V4.sh $Dir/Results/Final_Results_$z.tsv &
$ScriptPath/LowestCommonAncestor_V4.sh $Dir/Results/Final_Results_$z.tsv_AllResults.tsv &
wait
awk 'NR == 1; NR > 1 {print $0 | "sort -t',' -k3,3rn -k4,4rn"}' $Dir/Results/Final_Results_$z-lca.csv > $Dir/Results/Final_Results_$z-lca.sorted.csv  & #Normal sort, but it keeps the header at the top
awk 'NR == 1; NR > 1 {print $0 | "sort -t',' -k3,3rn -k4,4rn"}' $Dir/Results/Final_Results_${z}.tsv_AllResults-lca.csv > $Dir/Results/Final_Results_${z}_AllResults-lca.sorted.csv  & #Normal sort, but it keeps the header at the top
wait

### Generate wordcloud of species frequency and species tree
bsub -K -q C python2.7 $ScriptPath/CSV-to-WordCloudFormat.py $Dir/Results/Final_Results_$z-lca.sorted.csv $Dir/Results/$z.WordCloud.R $Dir/Results/$z.Wordcloud.pdf & #Create R script for wordcloud creation
bsub -K -q C python2.7 $ScriptPath/CSV-to-Tree.py $Dir/Results/Final_Results_$z-lca.sorted.csv $Dir/Results/$z.gv & #Create script for taxonomical tree creation
wait

### Gather all reads classified by BLAST
cat $Dir/Processing/ReadNames*.fsa > $Dir/Processing/All-Reads-From-BLAST_$z.fsa

### Print reads with predictions to file
bsub -K -q C python2.7 $ScriptPath/ReadNames-to-FASTA_V8.py $Dir/Results/Final_Results_$z.tsv_AllResults.tsv $Dir/Processing/All-Reads-From-BLAST_$z.fsa $Dir/${z}_v0.23_FungalReads.tsv_AllResults.fsa & 
bsub -K -q C python2.7 $ScriptPath/ReadNames-to-FASTA_V8.py $Dir/Results/Final_Results_$z.tsv $Dir/Processing/All-Reads-From-BLAST_$z.fsa $Dir/${z}_v0.23_FungalReads.fsa &

### Create text summary of pipeline
ClassifiedReads=$(wc -l $Dir/Processing/AllClassified_$z | awk '{print $1}')
RemovedReads="$((ReadsLeft - ClassifiedReads))"
echo "Number of reads removed by Kraken: " >> $PreDir/Run_Statistics.txt
echo $RemovedReads >> $PreDir/Run_Statistics.txt
echo "Number of reads remaining after Kraken: " >> $PreDir/Run_Statistics.txt
echo $ClassifiedReads >> $PreDir/Run_Statistics.txt
FungalReadsBLAST_FungiBacteria=$(wc -l $Dir/Results/Final_Results_$z.tsv | awk '{print $1}')
RemovedReadsBLAST_FungiBacteria="$((ClassifiedReads - FungalReadsBLAST_FungiBacteria))"
echo "Number of bacterial reads removed by BLAST: " >> $PreDir/Run_Statistics.txt
echo $RemovedReadsBLAST_FungiBacteria >> $PreDir/Run_Statistics.txt
echo "Number of reads predicted as fungal: " >> $PreDir/Run_Statistics.txt
echo $FungalReadsBLAST_FungiBacteria >> $PreDir/Run_Statistics.txt
FinishTime="Pipeline finished at $(date)"
echo -e $StartTime"\n"$FinishTime"\n" >> $PreDir/SummaryFile.txt

### Rename output files with fungal percentage as a prefix
FungalRead=$(bc -l <<< $FungalReadsBLAST_FungiBacteria/$ReadsIn)
FungalReadPercentage=$(python -c "print(float($FungalRead*100))")  #Get the percentage of the original reads that are predicted to be fungal
FungalReadPercentage2=$(printf "%0.2f\n" $FungalReadPercentage) #Shorten float to two decimal places
cp $Dir/Results/Final_Results_$z-lca.sorted.csv $Dir/${FungalReadPercentage2}_${z}.Results.csv
cp $Dir/Results/Final_Results_${z}_AllResults-lca.sorted.csv $Dir/${FungalReadPercentage2}_${z}.Results.IncludingFalsePositives.csv

### PDF creation
echo "Not enough fungal reads to create wordcloud" >> $Dir/Results/WordcloudError.txt
enscript -B -o $Dir/Results/WordcloudError.ps -f Times-Roman12 $Dir/Results/WordcloudError.txt
ps2pdf $Dir/Results/WordcloudError.ps $Dir/Results/$z.Wordcloud.pdf
dot -Tps $Dir/Results/$z.gv -o $Dir/Results/$z.TreeLineage.ps
gs -sDEVICE=pdfwrite -sOutputFile=$Dir/Results/$z.TreeLineage.pdf -dBATCH -dNOPAUSE $Dir/Results/$z.TreeLineage.ps
pdfcrop $Dir/Results/$z.TreeLineage.pdf $Dir/Results/Cropped.$z.TreeLineage.pdf
Rscript-3.3.3 $Dir/Results/$z.WordCloud.R
cat $FungTaxDir/PipelineSummary.txt $PreDir/SummaryFile.txt $PreDir/Run_Statistics.txt > $PreDir/FindFungi-TextFile.txt
enscript -B -o $PreDir/FindFungi-TextFile.ps -f Times-Roman12 $PreDir/FindFungi-TextFile.txt
ps2pdf $PreDir/FindFungi-TextFile.ps $PreDir/FindFungi-TextFile.pdf
enscript -B -o $Dir/Results/Final_Results_$z-lca.sorted.Fungi.ps -f Times-Roman12 $Dir/Results/Final_Results_$z-lca.sorted.csv
ps2pdf $Dir/Results/Final_Results_$z-lca.sorted.Fungi.ps $Dir/Results/Final_Results_$z-lca.sorted.Fungi.pdf
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFFitToPage -dPDFSETTINGS=/prepress -sOutputFile=$Dir/${FungalReadPercentage2}_${z}.Results.pdf $PreDir/FindFungi-TextFile.pdf $Dir/Results/$z.Wordcloud.pdf $Dir/Results/Cropped.$z.TreeLineage.pdf $Dir/Results/Final_Results_$z-lca.sorted.Fungi.pdf 
wait

### Clean up
mv *.std* $Dir/bsub_reports/
#rm -rf $PreDir  #Uncomment to remove all intermediate files. Keep in mind you will need to copy or move the results files to a safe directory prior to this step.

fi

echo "Finished at $(date)"

