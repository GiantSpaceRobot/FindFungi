### Download all fungal genomes from NCBI
### Gather representative/reference genomes
### Shuffle the genomes and split the FASTA files up
### Build Kraken databases

numberOfDBs=32



### Download all fungi from NCBI genbank
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes//genbank/fungi/* fungi/
find ./ -type f -exec chmod 644 {} \;  # Change mod rights of all NCBI files


### Gather representative genomes (we don't want all genomes as many assemblies are terrible)
python Genome_Grabber.py


### Shuffle sequences, DUST the genomes to mask repetitive regions
./Genome_DUSTer.sh
cat DUSTed_RepresentativeSpecies/* > All_DUSTed_RepresentativeSpecies.fsa # Add all dusted genomes to a single file
cat All_DUSTed_RepresentativeSpecies.fsa | perl seq-shuf.pl > All_DUSTed_RepresentativeSpecies_shuffled.fsa
pyfasta split -n $numberOfDBs All_DUSTed_RepresentativeSpecies_shuffled.fsa


### Loop over FASTA files in directory and create kraken databases out of each 
echo "Creating $numberOfDBs kraken databases..."
for i in $(seq 1 $numberOfDBs) 
    do
    mkdir -p chunks_$i # Create directories if they do not exist already
done       
     
for i in $(seq 1 $numberOfDBs)
    do
	### Copy taxonomy files from pre-existing database or download from NCBI
    cp -r /home/paul/Project-work/Databases/Kraken-Databases/Simulated-Metagenome-Original/taxonomy/ /home/paul/Project-work/Databases/Kraken-Databases/Kraken_32DB/chunks_$i & 
    ### Point to FASTA:
	kraken-build --threads 30 --add-to-library /home/paul/Project-work/Databases/Fungal-Genomes/NCBI-Fungi/Split-FASTA-EqualChunks/Chunks$i/Chunks$i.fsa --db chunks_$i &
done
wait
 
for i in $(seq 1 $numberOfDBs)
    do
    kraken-build --build --threads 30 --db chunks_$i &
done
wait
 
echo "Done"


