#!/usr/bin/bash

# Download and create BLAST database for Bant_MLVA31_analyzer.sh script

path=$(realpath ./)
DIR_PATH=$(realpath ./)
patterndir=$DIR_PATH/Supporting_files

cd $patterndir
mkdir $patterndir/Blast_dbase
cd $patterndir/Blast_dbase
wget -c https://zenodo.org/records/18399429/files/Bcereus_genomes_v3.fasta.gz
gzip -d *.gz
makeblastdb -in Bcereus_genomes_v3.fasta -dbtype nucl -parse_seqids -out Bcereus_genomes -title "Bacillus cereus group genomes"
