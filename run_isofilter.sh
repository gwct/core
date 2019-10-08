#!/bin/bash

############################################################
# For rodent comparative project, 10.19
# After retrieving sequences from Ensembl with 
# ensembl_get_fasta.sh, use this to retrieve the longest
# isoform of each gene.
############################################################

display_usage() {
# Usage message for -h
    echo
	echo "Usage:"
    echo
    echo "Read file directly:"
    echo "run_isofilter.sh -f [file] -o [both] -h"
    echo
    echo "Read from STDIN:"
    echo "cat [file with species names, one per line] | run_isofilter.sh -o [both] -h"
    echo
    echo "-f: An input file with one Ensembl species name per line. Can also be read from STDIN."
    echo "-o: The type of sequence to process: cds (coding sequences), pep (peptides), or both (default)."
    echo "-h: Display this help message."
    echo
    echo "====================================="
    exit 0
}

##########

invalid_opt() {
# Invalid input option message
    echo "Error 1: Invalid option: $1"
    echo "====================================="
    exit 1
}

##########
# Invalid output option message
invalit_out(){
    echo "Error 2: Invalid output (-o) option: $1"
    echo "Error 2: Output (-o) must be one of: both, pep, cds"
    echo "====================================="
    exit 1   
}

############################################################

echo
echo "====================================="
echo "Filter short isoforms from Ensembl sequences."
input="/dev/stdin"
output="both"
procs=1
# Defaults

while getopts ":f:o:p:h" arg; do
    case $arg in
        h) display_usage;;
        f) input=$OPTARG;;
        o) output=$OPTARG;;
        p) procs=$OPTARG;;
        \?) invalid_opt $OPTARG ;;
    esac
done
# Parse input options

if [ $output != "both" ] && [ $output != "pep" ] && [ $output != "cds" ]; then
    invalid_out $output
fi
# Check output option

echo "Reading species names from: " $input
echo "====================================="
echo
# Print info

while read spec; do
    if [[ $spec == \#* ]]; then
        continue
    fi
    echo $spec
    mkdir $spec
    #cd $spec
    if [ $output == "cds" ] || [ $output == "both" ]; then
        isofilter.py -i $spec/*.cds.all.fa.gz -t ens -o $spec/$spec-filtered-cds.fa
    fi

    if [ $output == "pep" ] || [ $output == "both" ]; then
        isofilter.py -i $spec/*.pep.all.fa.gz -t ens -o $spec/$spec-filtered-pep.fa
    fi
done < "$input"
# Run rsync

echo
echo "Done!"
echo "====================================="