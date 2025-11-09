#!/usr/bin/bash

# Download and create BLAST database for Bant_MLVA31_analyzer.sh script

path=$(realpath ./)
DIR_PATH=$(realpath ./)
patterndir=$DIR_PATH/Supporting_files

mkdir $DIR_PATH/Output_directory
cd $patterndir
mkdir $patterndir/Blast_dbase
cd $patterndir/Blast_dbase
wget -c https://zenodo.org/records/17078231/files/Bcereus_genomes_v2.fasta.gz
gzip -d *.gz
makeblastdb -in Bcereus_genomes_v2.fasta -dbtype nucl -parse_seqids -out Bcereus_genomes -title "Bacillus cereus group genomes"
