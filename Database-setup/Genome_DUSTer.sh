# DUST genomes to mask repetitive regions
 
for d in RepresentativeSpecies/*; do
    File=$(basename $d)
    mkdir -p DUSTed_RepresentativeSpecies
	dustmasker -outfmt fasta -in $d -out DUSTed_RepresentativeSpecies/$File
    echo "DUSTed $File"
    done

