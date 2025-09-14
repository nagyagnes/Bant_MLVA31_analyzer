#!/usr/bin/bash

# Testing of the installation of dependencies and running of Bant_MLVA31_analyzer.sh script 

path=$(readlink -f ${BASH_SOURCE:-$0})
DIR_PATH=$(dirname $path)

./Bant_MLVA31_analyzer.sh -s Test_1,Test_2 -o $DIR_PATH/Demo_files -m Y -i N
mv {$DIR_PATH/Demo_files/Test_1/Test_1.fastq,$DIR_PATH/Demo_files/Test_2/Test_2.fastq} $DIR_PATH/Demo_files
awk -F, 'NR==FNR { keys[$1]=1; next } !($1 in keys) { print $0 }' $DIR_PATH/Demo_files/Demo_results/Test_1/Test_1_final_results.csv $DIR_PATH/Demo_files/Test_1/Results/Test_1_final_results.csv > $DIR_PATH/Demo_files/Comparison_Test_1_results.txt
awk -F, 'NR==FNR { keys[$1]=1; next } !($1 in keys) { print $0 }' $DIR_PATH/Demo_files/Demo_results/Test_2/Test_2_final_results.csv $DIR_PATH/Demo_files/Test_2/Results/Test_2_final_results.csv > $DIR_PATH/Demo_files/Comparison_Test_2_results.txt
cat $DIR_PATH/Demo_files/Comparison_Test_1_results.txt $DIR_PATH/Demo_files/Comparison_Test_2_results.txt > $DIR_PATH/Demo_files/Comparison_test_results.txt
if [ -s $DIR_PATH/Demo_files/Comparison_test_results.txt ]
then
     echo "Test failed: differences found. Check dependencies."
else
     echo "Test successful."
     rm $DIR_PATH/Demo_files/Comparison_test_results.txt
fi
