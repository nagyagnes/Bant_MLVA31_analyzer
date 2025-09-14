#!/usr/bin/bash

# Author: Agnes Nagy (bio.lab.hu@gmail.com)
# version: 2025-08-26
#
# Copyright (c) 2025 Ágnes Nagy (bio.lab.hu@gmail.com)
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
############# SETUP OF OPTIONS, VARIABLES, INPUT FILES ###################
##########################################################################

usage="Bant_MLVA31_analyzer -- performs 31-loci MLVA typing for Bacillus anthracis directly from ONT reads from environmental samples, and finds the closest relative strains\n
Version 2025-08-26\n
\n
Options: [-h] [-s -p -d -o -m -i]\n
where:\n
-h=Prints this help text\n
-s=Set a name or number for sample which is the name of the input .fasq file [Required]\n
-p=Set the directory for supporting files [Default: Supporting_files]\n
-d=Set the directory for scripts [Default: Script_directory]\n
-o=Set the output directory which contains input sample.fastq file and will contain output directories for intermediate and result files [Default: Output_directory]\n
-m=Determine the closest neighbour strain from global Bacillus anthracis database [Required, Y=yes, N=no]\n
-i=Keep intermediate files and directories created during analysis [Required, Y=yes, N=no]"

now="$(date +'%Y-%m-%d %T')"
version="Bant_MLVA31_analyzer Version 2025-08-26"
path=$(readlink -f ${BASH_SOURCE:-$0})
DIR_PATH=$(dirname $path)
patterndir=$DIR_PATH/Supporting_files
scriptdir=$DIR_PATH/Script_directory
OUTPUT=$DIR_PATH/Output_directory
while getopts 's:p:d:o:m:i:h' options; do
    case "${options}" in
        h)  echo -e $usage && exit       ;;
        s)  set -f
            IFS=,
            strain=($OPTARG)	          ;;
        p)  patterndir="${OPTARG}"       ;;
        d)  scriptdir="${OPTARG}"        ;;
        o)  OUTPUT="${OPTARG}"	  ;;
        m)  set -f
            closest_neighbour=($OPTARG)  ;;
        i)  set -f
            intermediate_files=($OPTARG) ;;
   *)	echo -e $usage		          ;;
esac
done
shift $((OPTIND-1))

##########################################################################
################# CHECKING INPUT FILES AND DIRECTORIES ###################
##########################################################################
((

if [ -z $strain ]; then
	echo -e "\nError: No input sample was set !!\n"
	echo -e $usage
	exit;
fi

if [ -z $closest_neighbour ]; then
	echo -e "\nError: No option for closest neigbour determination was set !!\n"
	echo -e $usage
	exit;
fi

case "$closest_neighbour" in
		Y) echo "Closest neighbour determination is selected"                                     ;;
		N) echo "Closest neighbour determination is not selected"	                            ;;
		*) echo "Invalid selection for closest neighbour determination" && echo -e $usage && exit ;;
esac

if [ -z $intermediate_files ]; then
	echo -e "\nError: No option for keeping or deleting of intermediate files was set !!\n"
	echo -e $usage
	exit;
fi

case "$intermediate_files" in
		Y) echo "Keeping of intermediate files is selected"                                          ;;
		N) echo "Keeping of intermediate files is not selected. Intermediate files will be deleted!" ;;
		*) echo "Invalid selection for keeping of intermediate files" && echo -e $usage && exit      ;;
esac
	
echo -e "\nInput samples:"
for i in ${!strain[@]}; do echo "${strain[$i]}"; done

if [[ ! -d "$OUTPUT" ]]; then
	echo -e "\nError: Output directory does not exist !!\n"
	echo -e $usage
	exit;
fi

echo -e "\nOutput directory: "$OUTPUT""

if [[ ! -d "$scriptdir" ]]; then
	echo -e "\nError: Script directory does not exist !!\n"
	echo -e $usage
	exit;
fi

echo -e "\nScript directory: "$scriptdir""

cd $scriptdir
if [[ ! -f "Compute_Sample_Length_v3.R" ]] && [[ ! -f "MLVA_read_length_distribution.R" ]] && [[ ! -f "MLVA_distmatrix_v3_modMatchCoeffwoNA_closestrelative.R" ]] && [[ ! -f "Calculation_template.csv" ]] ; then
	echo -e "\nError: Supporting scripts do not exist !!\n"
	exit;
fi

cd $DIR_PATH
if [[ ! -d "$patterndir" ]]; then
	echo -e "\nError: Supporting files directory does not exist !!\n"
	echo -e $usage
	exit;
fi

echo -e "\nSupporting files directory: "$patterndir""

for i in ${!strain[@]}; do
	cd $OUTPUT
	if [[ ! -f "$OUTPUT/${strain[$i]}.fastq" ]]; then
		echo -e "\nError: ${strain[$i]}.fastq does not exist !!"
		exit;
	fi
done

ref=$patterndir/Reference/BaAmesAncestor.fasta
region=$patterndir/BaAmes_MLVA_regions.txt
blastdbase=$patterndir/Blast_dbase/Bcereus_genomes
repeats=("bams01" "bams03" "bams05" "bams13" "bams15" "bams21" "bams22" "bams23" "bams24" "bams25" "bams28" "bams30" "bams31" "bams34" "bams44" "bams51" "bams53" "CG3" "pXO1" "pXO2" "vntr12" "vntr16" "vntr17" "vntr19" "vntr23" "vntr35" "vrrA" "vrrB1" "vrrB2" "vrrC1" "vrrC2")
repeats1=("CG3" "pXO1" "pXO2" "vntr12" "vntr19" "vntr23" "vntr35")

##########################################################################
##################### START OF ANALYSIS###################################
##########################################################################
SECONDS=0

echo -e "\n### $version\n### Analysis started at $now\n"

### Set up a loop for multiple samples (.fastq files) to be processed ####
for i in ${!strain[@]}; do
	cd $OUTPUT
	## Create directories for input, intermediate and result files ##
	mkdir ${strain[$i]}
	workdir=$OUTPUT/${strain[$i]}
	mv ${strain[$i]}.fastq $workdir
	mkdir -p $workdir/{Map,Repeats,Consensus,Close_relatives,Results}
	mkdir -p $workdir/Consensus/{Reads_blast,Reads_sorted}
	echo "$workdir/Repeats" > $scriptdir/workdir2.txt
	echo "$workdir/Close_relatives" > $scriptdir/workdir3.txt
	
	##########################################################################
        ############## I. SEPARATION OF READS BY REPEAT REGIONS ##################
        ##########################################################################
	## Map reads to Bacillus anthracis reference genome sequence ##
	echo -e "\n## ${strain[$i]} map to reference started. ##\n"
	cd $workdir/Map
	minimap2 --MD -L -t 10 -ax map-ont $ref $workdir/${strain[$i]}.fastq -o ${strain[$i]}.sam
	samtools view -S -b ${strain[$i]}.sam -o ${strain[$i]}.bam
	samtools sort -o ${strain[$i]}.sorted.bam ${strain[$i]}.bam
	samtools index ${strain[$i]}.sorted.bam
	echo -e "\n## ${strain[$i]} map to reference completed. ##\n"
	
	## Set up a loop for processing all 31 repeats per sample ##
	for k in ${!repeats[@]}; do
		echo -e "## ${strain[$i]} ${repeats[$k]} analysis started. ##\n"
		cd $workdir/Repeats
		## Split the mapped reads by the repeat regions of the reference genome ##
		IFS=$'\n' read -d '' -r -a position < $region
		samtools view -bo ${repeats[$k]}.bam -X $workdir/Map/${strain[$i]}.sorted.bam $workdir/Map/${strain[$i]}.sorted.bam.bai ${position[k]}
		## Make .fastq files from .bam ##
		samtools bam2fq ${repeats[$k]}.bam > ${repeats[$k]}.fastq
		echo "${repeats[$k]}" $(expr $(cat ${repeats[$k]}.fastq | wc -l) / 4) > ${repeats[$k]}_fastq_readnr.txt
		
		###############################################################################################################
		############## II. READ FILTRATION WITH ALIGNING TO COLLECTION OF REPEAT REFERENCE SEQUENCES ##################
		###############################################################################################################		
		## Map splitted reads to collection of repeat reference sequences ##
		cd $patterndir
		if ! test -f "${repeats[$k]}_ref.mmi"; then
			minimap2 -d ${repeats[$k]}_ref.mmi ${repeats[$k]}_ref.fa
			samtools faidx ${repeats[$k]}_ref.fa
		fi
		cd $workdir/Repeats	
		minimap2 --MD -L -t 10 -ax map-ont $patterndir/${repeats[$k]}_ref.mmi ${repeats[$k]}.fastq > ${repeats[$k]}_ref.sam
		samtools view -bS -T $patterndir/${repeats[$k]}_ref.fa ${repeats[$k]}_ref.sam > ${repeats[$k]}_ref.bam
		samtools sort ${repeats[$k]}_ref.bam -o ${repeats[$k]}_ref.sorted.bam
		samtools index ${repeats[$k]}_ref.sorted.bam
		## Select reference sequence with maximum mean depth of coverage ##
		samtools coverage ${repeats[$k]}_ref.sorted.bam -o ${repeats[$k]}_ref.csv
		gawk 'NR>1 {if(max<$7){max=$7;line1=$1; line2=$3; line3=$4}}END{print line1, line2, line3;}' ${repeats[$k]}_ref.csv > ${repeats[$k]}_ref_unit.csv
		unitnr=$(awk -F' ' '{print $1}' ${repeats[$k]}_ref_unit.csv)
		samtools view -bo ${repeats[$k]}_ref_unit.bam -X ${repeats[$k]}_ref.sorted.bam ${repeats[$k]}_ref.sorted.bam.bai $unitnr
		## Export reads mapped to selected reference sequence ##
		samtools bam2fq ${repeats[$k]}_ref_unit.bam > ${repeats[$k]}_ref_unit.fastq
		
		###############################################################################################################
		############## III. IN SILICO MLVA31 PCR OF READS, FILTRATION BY AMPLICON LENGTH DISTRIBUTION #################
		##############               AND DETERMINATION OF THE AMPLICON LENGTH                         #################
		###############################################################################################################	
		## Search for amplicons in reads with primersearch ##
		primersearch -seqall ${repeats[$k]}_ref_unit.fastq -infile $patterndir/Primer_files/${repeats[$k]}_primer -mismatchpercent 10 -outfile ${repeats[$k]}_ref_unit.primers.primersearch
		## Export the length value of the most frequent amplicons with the same length  ##
		perl -i -pe 's/[ \t]+//g' ${repeats[$k]}_ref_unit.primers.primersearch
		perl -i -n -e 'print if /\S/' ${repeats[$k]}_ref_unit.primers.primersearch
		grep -oP '(?<=Amplimerlength:).*(?=bp)' ${repeats[$k]}_ref_unit.primers.primersearch > ${repeats[$k]}_ref_unit.sizes.txt
		## Determine the number of reads containing amplicon with exported length value ##
		readnr=$(wc -l < ${repeats[$k]}_ref_unit.sizes.txt)
		echo -e "\n${repeats[$k]} number of selected amplicons: "$readnr"\n"
		## Set up an if statement according to determined number of reads ##
		if (($readnr > 4)); then
		
		###############################################################################################################
		########################## III/1. READ FILTRATION BY AMPLICON LENGTH DISTRIBUTION #############################
		###############################################################################################################
			## Compute the descriptive statistics of the length distribution of the amplimers with R ##
			echo "${repeats[$k]}_ref_unit.sizes.txt" > workfile.txt
			cd $scriptdir
			Rscript MLVA_read_length_distribution.R
			## Export list of reads containing amplicon length between minimum and maximum value ##
			cd $workdir/Repeats
			rm workfile.txt
			awk '{if (NR!=1) {print}}' ${repeats[$k]}_ref_unit.primers.primersearch | awk 'ORS=(NR%5==0)?"\n":"\t"' | awk '{print $2, $5}' > ${repeats[$k]}_ref_unit.amplimers.txt
			min=$(awk -F';' 'NR>1 {print $1}' ${repeats[$k]}_ref_unit.stat.txt)
			max=$(awk -F';' 'NR>1 {print $2}' ${repeats[$k]}_ref_unit.stat.txt)
			gawk -F 'Sequence:|Amplimerlength:|bp' '{print $2, $3}' ${repeats[$k]}_ref_unit.amplimers.txt | sort -k2 | awk -v a="$min" -v b="$max" '$2 >= a && $2 <= b {print $1}' > ${repeats[$k]}_ref_unit.lst
			## Determine the number of reads in the list ##
			readnr2=$(wc -l < ${repeats[$k]}_ref_unit.lst)
			echo -e "\nRead number after amplicon length filtering:" "$readnr2"
			echo "${repeats[$k]}" "$readnr2" > ${repeats[$k]}_lst_readnr.csv
			## Set up an if statement according to determined number of reads ##
			if (($readnr2 > 20)); then
			
			###############################################################################################################
		        ############## III/1/1. CONSENSUS SEQUENCE CALLING, CHECKING SIMLIRATITY TO BACILLUS ANTHRACIS ################
		        ###############################################################################################################
				## Export reads in read list to .fastq file ##
				seqtk subseq ${repeats[$k]}_ref_unit.fastq ${repeats[$k]}_ref_unit.lst > ${repeats[$k]}_ref_unit_reads.fastq
				## Generate consensus sequence(s) with amplicon_sorter ##
				echo -e "\n## ${repeats[$k]} consensus sequence calling started. ##\n"
				cd $scriptdir
				python3 amplicon_sorter.py -i $workdir/Repeats/${repeats[$k]}_ref_unit_reads.fastq -o $workdir/Consensus/ -np 4 -min 190 -max 3000 -maxr 1000
				find . -name "*.pickle" -exec rm {} \;
				find . -name "*.tmp" -exec rm {} \;
				cd $workdir/Consensus
				if [ -d $workdir/Consensus/${repeats[$k]}_ref_unit_reads ]; then
					cd $workdir/Consensus/${repeats[$k]}_ref_unit_reads
					mv ${repeats[$k]}_ref_unit_reads_consensussequences.fasta $workdir/Consensus/
				else
					echo -e "\n## ${repeats[$k]} consensus sequence calling failed. ##\n"
				fi
				## Trim the consensus sequence(s) with cutadapt retaining primers on both ends of sequence ##
				cd $workdir/Consensus
				primers=$(awk -F' ' '{print $1}' $patterndir/Primer_files/${repeats[$k]}_primercut.txt)
				if test -f "${repeats[$k]}_ref_unit_reads_consensussequences.fasta"; then
					cutadapt -g $primers --action=retain --overlap=10 --revcomp -o ${repeats[$k]}_cons_final.fasta ${repeats[$k]}_ref_unit_reads_consensussequences.fasta -j 8
				else
					echo "NA" > ${repeats[$k]}_cons_final.fasta
				fi
				## Align consensus sequence(s) to NCBI Bacillus cereus genome database with offline blast ##
				blastn -db $blastdbase -query ${repeats[$k]}_cons_final.fasta -task megablast -dust no -outfmt "7 delim=, qacc sacc bitscore score evalue pident nident qstart qend sstart send" -max_target_seqs 10 -subject_besthit -num_threads 8 > ${repeats[$k]}_cons_blast.csv
				## Select the name of the subject sequence with the highest percentage of identical matches and bitscore from blast results ##
				sed -i '/#/d' ${repeats[$k]}_cons_blast.csv
				sort -t, -k 6nr -k 3nr -k 2nr ${repeats[$k]}_cons_blast.csv | awk -F, 'NR == 1 {print $1, $2}' > ${repeats[$k]}_cons_seq_blast.csv
				echo "${repeats[$k]}" "Consensus_Bant" > ${repeats[$k]}_comment.csv
				cd $workdir/Consensus
				species=$(awk -F' ' '{print substr ($2,1,6)}' ${repeats[$k]}_cons_seq_blast.csv)
				## Set up a case statement based on results of blast ##
				case "$species" in
				Banthr) ## Determine the exact lenght of repeat region by the length of generated consensus sequence ##
					cd $workdir/Consensus
					sed 's/>//g' ${repeats[$k]}_cons_final.fasta > ${repeats[$k]}_cons_final.txt
					cons=$(awk -F' ' '{print $1}' ${repeats[$k]}_cons_seq_blast.csv)
					sed 's/).*/)/' ${repeats[$k]}_cons_final.txt | awk 'ORS=NR%2?" ":"\n"' | awk -v c="$cons" '$1 == c {print $2}' > ${repeats[$k]}_cons_final_seq.txt
					awk -v name="${repeats[$k]}" '{print name, length}' ${repeats[$k]}_cons_final_seq.txt > $workdir/Repeats/${repeats[$k]}_ref_unit.size_final.txt
					echo -e ">${strain[$i]}_${repeats[$k]}\n$(cat ${repeats[$k]}_cons_final_seq.txt)" > $workdir/Results/${strain[$i]}_${repeats[$k]}_cons.fasta
					;; ## 1. OUTPUT OF III. step ##
				*)
				###############################################################################################################
				################### III/1/1/1. SIMILARITY CHECKING OF FILTERED READS TO BACILLUS ANTHRACIS ####################
				###############################################################################################################	
					cd $workdir/Consensus
					[ -f ${repeats[$k]}_cons_final.fasta ] && rm ${repeats[$k]}_cons_final.fasta
					[ -f ${repeats[$k]}_cons_seq_blast.csv ] && rm ${repeats[$k]}_cons_seq_blast.csv
					[ -f ${repeats[$k]}_comment.csv ] && rm ${repeats[$k]}_comment.csv
					cd $workdir/Repeats
					## Export reads listed in III/1. step to .fastq file ##
					seqtk subseq ${repeats[$k]}_ref_unit.fastq ${repeats[$k]}_ref_unit.lst > ${repeats[$k]}_ref_unit_reads.fastq
					## Convert fastq files to fasta ##
					seqtk seq -A ${repeats[$k]}_ref_unit_reads.fastq > ${repeats[$k]}_ref_unit_reads.fasta
					cd $workdir/Consensus
					echo "${repeats[$k]}" "No_cons_Bant" > ${repeats[$k]}_comment.csv
					## Trim the reads retaining primers on both ends of sequence ##
					primers=$(awk -F' ' '{print $1}' $patterndir/Primer_files/${repeats[$k]}_primercut.txt)
					cutadapt -g $primers --action=retain --overlap=10 --revcomp -o ${repeats[$k]}_ref_unit_reads_cutted.fasta $workdir/Repeats/${repeats[$k]}_ref_unit_reads.fasta -j 8
					## Align read sequences to NCBI Bacillus cereus genome database ##
					blastn -db $blastdbase -query ${repeats[$k]}_ref_unit_reads_cutted.fasta -task megablast -dust no -outfmt "7 delim=, qacc sacc bitscore pident" -max_target_seqs 10 -subject_besthit -num_threads 8 > $workdir/Consensus/Reads_blast/${repeats[$k]}_blast_reads_cutted.csv
					## Select the name of the subject with the highest percentage of identical matches and bitscore from blast results of each read ##
					cd $workdir/Consensus/Reads_blast
					awk -F'\t' '$0 !~ "^(#)"' ${repeats[$k]}_blast_reads_cutted.csv > ${repeats[$k]}_blast_reads_cuttedf.csv
					awk -F ',' -v name="${repeats[$k]}" '{print >> ($1"."name".csv")}' ${repeats[$k]}_blast_reads_cuttedf.csv
					mv ${repeats[$k]}_blast_reads_cutted.csv ${repeats[$k]}_blast_reads_cuttedf.csv $workdir/Consensus
					find . -name "*.${repeats[$k]}.csv" -print0 | while read -d $'\0' file
					do
						sort -t, -k 4nr -k 3nr -k 2nr "$file" | head -1 > "$workdir/Consensus/Reads_sorted/$file";  
					done
					cd $workdir/Consensus/Reads_sorted
					find . -name "*.${repeats[$k]}.csv" -exec cat {} > $workdir/Consensus/${repeats[$k]}_reads_blast_result.csv \;
					cd $workdir/Consensus
					## Set up an if statement based on results of reads blast ##
					Banthread=$(wc -l < ${repeats[$k]}_reads_blast_result.csv)
					Allread=$(awk 'END { print NR - 1 }' RS='Banthr' ${repeats[$k]}_reads_blast_result.csv)
					Readpercentage=$(gawk -v a="$Banthread" -v b="$Allread" 'BEGIN { x = a; y = b; z = 100; print ((y/x)*z)}')
					R=${Readpercentage%%.*}
					if (($R > 90)); then
						## Determine the exact lenght of repeat using the value of median from descriptive statistics computed in III/1. step ##
						cd $workdir/Repeats
						awk -F';' -v name="${repeats[$k]}" 'NR>1 {print name, $3}' ${repeats[$k]}_ref_unit.stat.txt > ${repeats[$k]}_ref_unit.size_final.txt
						## 2. OUTPUT OF III. step ##
					else
						## The length of repeat region can not be determined ##
						cd $workdir/Repeats
						echo "${repeats[$k]}" "NA" > ${repeats[$k]}_ref_unit.size_final.txt
						cd $workdir/Consensus
						rm ${repeats[$k]}_reads_blast_result.csv
						echo "NA" > ${repeats[$k]}_reads_blast_result.csv
						echo "${repeats[$k]}" "No_cons_NoBant" > ${repeats[$k]}_comment.csv
						## 3. OUTPUT OF III. step ##
					fi
				esac
			else
			###############################################################################################################
			#################### III/1/2. SIMILARITY CHECKING OF FILTERED READS TO BACILLUS ANTHRACIS #####################
			###############################################################################################################	
				## Export reads listed in III/1. step to .fastq file ##
				cd $workdir/Repeats
				seqtk subseq ${repeats[$k]}_ref_unit.fastq ${repeats[$k]}_ref_unit.lst > ${repeats[$k]}_ref_unit_reads.fastq
				## Convert fastq files to fasta ##
				seqtk seq -A ${repeats[$k]}_ref_unit_reads.fastq > ${repeats[$k]}_ref_unit_reads.fasta
				cd $workdir/Consensus
				echo "${repeats[$k]}" "No_cons_LowReadBant" > ${repeats[$k]}_comment.csv
				echo "${repeats[$k]}" "$readnr" > $workdir/Repeats/${repeats[$k]}_lst_readnr.csv
				## Trim the reads retaining primers on both ends of sequence ##
				primers=$(awk -F' ' '{print $1}' $patterndir/Primer_files/${repeats[$k]}_primercut.txt)
				cutadapt -g $primers --action=retain --overlap=10 --revcomp -o ${repeats[$k]}_ref_unit_reads_cutted.fasta $workdir/Repeats/${repeats[$k]}_ref_unit_reads.fasta -j 8
				## Align read sequences to NCBI Bacillus cereus genome database ##
				blastn -db $blastdbase -query ${repeats[$k]}_ref_unit_reads_cutted.fasta -task megablast -dust no -outfmt "7 delim=, qacc sacc bitscore pident" -max_target_seqs 10 -subject_besthit -num_threads 8 > $workdir/Consensus/Reads_blast/${repeats[$k]}_blast_reads_cutted.csv
				## Select the name of the subject with the highest percentage of identical matches and bitscore from blast results of each read ##
				cd $workdir/Consensus/Reads_blast
				awk -F'\t' '$0 !~ "^(#)"' ${repeats[$k]}_blast_reads_cutted.csv > ${repeats[$k]}_blast_reads_cuttedf.csv
				awk -F ',' -v name="${repeats[$k]}" '{print >> ($1"."name".csv")}' ${repeats[$k]}_blast_reads_cuttedf.csv
				mv ${repeats[$k]}_blast_reads_cutted.csv ${repeats[$k]}_blast_reads_cuttedf.csv $workdir/Consensus
				find . -name "*.${repeats[$k]}.csv" -print0 | while read -d $'\0' file
				do
					sort -t, -k 4nr -k 3nr -k 2nr "$file" | head -1 > "$workdir/Consensus/Reads_sorted/$file";  
				done
				cd $workdir/Consensus/Reads_sorted
				find . -name "*.${repeats[$k]}.csv" -exec cat {} > $workdir/Consensus/${repeats[$k]}_reads_blast_result.csv \;
				## Set up an if statement based on results of reads blast ##
				cd $workdir/Consensus
				Banthread=$(wc -l < ${repeats[$k]}_reads_blast_result.csv)
				Allread=$(awk 'END { print NR - 1 }' RS='Banthr' ${repeats[$k]}_reads_blast_result.csv)
				Readpercentage=$(gawk -v a="$Banthread" -v b="$Allread" 'BEGIN { x = a; y = b; z = 100; print ((y/x)*z)}')
				R=${Readpercentage%%.*}
				if (($R > 50)); then
				        ## Determine the exact lenght of repeat using the value of median from descriptive statistics computed in III/1. step ##
					cd $workdir/Repeats
					awk -F';' -v name="${repeats[$k]}" 'NR>1 {print name, $3}' ${repeats[$k]}_ref_unit.stat.txt > ${repeats[$k]}_ref_unit.size_final.txt
					## 4. OUTPUT OF III. step ##
				else
				        ## The length of repeat region can not be determined ##	
					cd $workdir/Repeats
					echo "${repeats[$k]}" "NA" > ${repeats[$k]}_ref_unit.size_final.txt
					cd $workdir/Consensus
					echo "${repeats[$k]}" "No_cons_LowReadNoBant" > ${repeats[$k]}_comment.csv
					## 5. OUTPUT OF III. step ##
				fi
			fi
		else ## The length of repeat region can not be determined ##
			cd $workdir/Repeats
			echo "${repeats[$k]}" "NA" > ${repeats[$k]}_ref_unit.size_final.txt
			echo "${repeats[$k]}" "NA" > ${repeats[$k]}_lst_readnr.csv
			cd $workdir/Consensus
			echo "NA" > ${repeats[$k]}_cons_seq_blast.csv
			echo "${repeats[$k]}" "No_Read" > ${repeats[$k]}_comment.csv
			## 6. OUTPUT OF III. step ##
		fi
		echo -e "\n## ${strain[$i]} ${repeats[$k]} analysis completed. ##\n"
	done
	echo -e "### ${strain[$i]}  sample all repeats completed. ###\n"
	#####################################################################################
        ############## IV. REPEATS > 7 BASES: CORRECTION OF REPEAT LENGTHS ##################
        #####################################################################################
	## If the blast result is Bacillus anthracis correct the length of repeat region for repeats shorter than 7 basepairs using the selected reference sequence in II. step ##
	echo -e "### ${strain[$i]} Short repeat length correction started. ###\n"
	for m in ${!repeats1[@]}; do
		cd $workdir/Consensus
		if test -f "${repeats1[$m]}_cons_seq_blast.csv"; then
			species=$(awk -F' ' '{print substr ($2,1,6)}' ${repeats1[$m]}_cons_seq_blast.csv)
			case "$species" in
			Banthr)
				cd $workdir/Repeats
				rm ${repeats1[$m]}_ref_unit.size_final.txt
				awk -F' ' -v name="${repeats1[$m]}" '{print name, $2}' ${repeats1[$m]}_ref_unit.csv > ${repeats1[$m]}_ref_unit.size_final.txt
				;;
			*)
				cd $workdir/Repeats
				rm ${repeats1[$m]}_ref_unit.size_final.txt
				echo "${repeats1[$m]}" "NA" > ${repeats1[$m]}_ref_unit.size_final.txt
			esac
		fi
		if test -f "${repeats1[$m]}_reads_blast_result.csv"; then
			Banthread=$(wc -l < ${repeats1[$m]}_reads_blast_result.csv)
			Allread=$(awk 'END { print NR - 1 }' RS='Banthr' ${repeats1[$m]}_reads_blast_result.csv)
			Readpercentage=$(gawk -v a="$Banthread" -v b="$Allread" 'BEGIN { x = a; y = b; z = 100; print ((y/x)*z)}')
			R=${Readpercentage%%.*}
			if (($R > 50)); then
				cd $workdir/Repeats
				rm ${repeats1[$m]}_ref_unit.size_final.txt
				awk -F' ' -v name="${repeats1[$m]}" '{print name, $2}' ${repeats1[$m]}_ref_unit.csv > ${repeats1[$m]}_ref_unit.size_final.txt
			else
				cd $workdir/Repeats
				rm ${repeats1[$m]}_ref_unit.size_final.txt
				echo "${repeats1[$m]}" "NA" > ${repeats1[$m]}_ref_unit.size_final.txt
			fi
		fi
		echo -e "# ${strain[$i]} ${repeats1[$m]} repeat length correction completed. #\n"
	done
	cd $OUTPUT
	echo -e "### ${strain[$i]} Short repeat length correction completed. ###\n"
	#####################################################################################
        ####################### V. CALCULATION OF FINAL RESULTS #############################
        #####################################################################################
        ## CONCATENATE OUTPUT FILES OF REPEATS ##
	cd $workdir/Repeats
	find . -name "*_ref_unit.size_final.txt" -exec cat {} > $workdir/${strain[$i]}.MLVA_sizes.txt \;
	find . -name "*_lst_readnr.csv" -exec cat {} > ${strain[$i]}.lst_readnr_sum.csv \;
	find . -name "*_fastq_readnr.txt" -exec cat {} > ${strain[$i]}_fastq_readnr_sum.csv \;
	cd $workdir/Consensus
	find . -name "*_comment.csv" -exec cat {} > $workdir/${strain[$i]}.comment.txt \;
	echo $workdir > $scriptdir/workdir.txt
	## COMPUTE REPEAT NUMBERS AND SAVE FINAL RESULT FILE ##
	echo -e "\n### ${strain[$i]} repeat number calculation started. ###\n"
	cd $scriptdir
	Rscript Compute_Sample_Length_v3.R
	rm workdir.txt workdir2.txt readme.txt
	cd $workdir
	find . -name "*.MLVA_unit.csv" -exec mv {} $workdir/Results \;
	mv ${strain[$i]}.MLVA_sizes.txt ${strain[$i]}.comment.txt $workdir/Repeats
	echo -e "### ${strain[$i]} repeat number calculation finished. ###\n"
	######################################################################################
	################# VI. Determination of closest relative strains ######################
	######################################################################################
	case "$closest_neighbour" in
		Y)
			cd $workdir/Results
			find . -type f -name ${strain[$i]}.MLVA_unit.csv -exec mv {} $workdir/Close_relatives \;
			cd $workdir/Close_relatives
			awk -F';' '{print $3}' ${strain[$i]}.MLVA_unit.csv > ${strain[$i]}.unitnr.csv
			sed -i "1s/.*/${strain[$i]}/" ${strain[$i]}.unitnr.csv
			sed -i -z 's/\n/;/g' ${strain[$i]}.unitnr.csv
			sed -i 's/.$//' ${strain[$i]}.unitnr.csv
			NAnr=$(awk '{print gsub(/\<NA\>/, "")}' ${strain[$i]}.unitnr.csv)
			echo -e "\n### ${strain[$i]} closest relative determination started. ###\n"
			if (($NAnr < 31)); then
				cat $patterndir/BaMLVA31_genotypes_R.csv ${strain[$i]}.unitnr.csv > ${strain[$i]}_BaMLVA31_genotypes_R.csv
				echo "${strain[$i]}" > samplename.txt
				echo "${strain[$i]}_BaMLVA31_genotypes_R.csv" > filename.txt
				cd $scriptdir
				Rscript MLVA_distmatrix_v3_modMatchCoeffwoNA_closestrelative.R
				rm workdir3.txt
				cd $workdir/Close_relatives
				rm filename.txt
				rm samplename.txt
				sed -i 's/"//g' ${strain[$i]}_close_rel.csv
				awk -i inplace -F';' 'NR>1 {print}' ${strain[$i]}_close_rel.csv
				sed -i '/^'"${strain[$i]}"'/d' ${strain[$i]}_close_rel.csv
				sort -k1 -o ${strain[$i]}_close_rel.csv{,}
				awk -F';' '{print $1}' ${strain[$i]}_close_rel.csv > ${strain[$i]}_close_rel.lst
				grep -Fwf ${strain[$i]}_close_rel.lst $patterndir/BaMLVA_genotypes_ref.csv > ${strain[$i]}_close_rel_data.csv
				sort -k1 -o ${strain[$i]}_close_rel_data.csv{,}
				paste -d';' ${strain[$i]}_close_rel_data.csv  <(awk -F';' '{print $2}' ${strain[$i]}_close_rel.csv) > ${strain[$i]}_close_rel_final.csv
				sort -t ';' -k35 -o ${strain[$i]}_close_rel_final.csv{,}
				header=$(awk -F';' '{print}' $patterndir/Headerfile.csv)
				sed -i '1i\'$header'\' ${strain[$i]}_close_rel_final.csv
				cat ${strain[$i]}.MLVA_unit.csv ${strain[$i]}_close_rel_final.csv > $workdir/Results/${strain[$i]}_final_results.csv
			else
				echo "Repeat numbers and closest relatives not determined" > ${strain[$i]}_close_rel_final.csv
				cat ${strain[$i]}.MLVA_unit.csv ${strain[$i]}_close_rel_final.csv > $workdir/Results/${strain[$i]}_final_results.csv
			fi
			;;
		N) 
			echo -e "\n### ${strain[$i]} closest relative determination is not selected. ###\n"
			cd $workdir/Results
			sed '$ a\Determination of closest relatives was not selected.' ${strain[$i]}.MLVA_unit.csv > ${strain[$i]}_final_results.csv
			rm ${strain[$i]}.MLVA_unit.csv
			;;
		esac
	case "$intermediate_files" in
		Y)
			echo -e "\n### ${strain[$i]} Keeping of intermediate directories and file is selected. Intermediate file in directories /Map, /Repeats, /Consensus, /Close_relatives saved. ###\n"
			;;
		N)
			cd $workdir
			rm -rf $workdir/Map $workdir/Repeats $workdir/Consensus $workdir/Close_relatives
			echo -e "\n### ${strain[$i]} Keeping of intermediate directories and file is not selected. Intermediate file in directories /Map, /Repeats, /Consensus, /Close_relatives deleted. ###\n"
			;;
		esac
	echo -e "\n### Analysis of ${strain[$i]} sample finished. ###\n"
done
echo -e "##### All Done. #####\n"
##########################################################################
####################### END OF ANALYSIS###################################
##########################################################################

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
) 2>&1) | tee -a Bant_MLVA31_analyzer.log
## END OF SCRIPT ##
