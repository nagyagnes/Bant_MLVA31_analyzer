#!/usr/bin/bash

# Download and create BLAST database for Bant_MLVA31_analyzer.sh script

path=$(realpath ./)
DIR_PATH=$(realpath ./)
patterndir=$DIR_PATH/Supporting_files
repeats=("bams01" "bams03" "bams05" "bams13" "bams15" "bams21" "bams22" "bams23" "bams24" "bams25" "bams28" "bams30" "bams31" "bams34" "bams44" "bams51" "bams53" "CG3" "pXO1" "pXO2" "vntr12" "vntr16" "vntr17" "vntr19" "vntr23" "vntr35" "vrrA" "vrrB1" "vrrB2" "vrrC1" "vrrC2")

cd $patterndir
mkdir $patterndir/Blast_dbase
cd $patterndir/Blast_dbase
wget -c https://zenodo.org/records/17078231/files/Bcereus_genomes_v2.fasta.gz
gzip -d *.gz
makeblastdb -in Bcereus_genomes_v2.fasta -dbtype nucl -parse_seqids -out Bcereus_genomes -title "Bacillus cereus group genomes"
cd $patterndir
for k in ${!repeats[@]}; do
	minimap2 -d ${repeats[$k]}_ref.mmi ${repeats[$k]}_ref.fa
	samtools faidx ${repeats[$k]}_ref.fa
done
