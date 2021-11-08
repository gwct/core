#!/bin/bash

############################################################
# Retrieves fasta files from Ensembl
############################################################

display_usage() {
# Usage message for -h
    echo
    echo "Usage:"
    echo
    echo "Read file directly:"
    echo "ensembl_get_orths.sh -r [98] -f [file] -o [both] -h"
    echo
    echo "Read from STDIN:"
    echo "cat [file with species names, one per line] | ensembl_get_orths.sh -r [98] -o [both] -h"
    echo
    echo "-r: Ensembl release. Default: 98"
    echo "-f: An input file with one Ensembl species name per line. Can also be read from STDIN. Lines starting with # are skipped."
    echo "-o: The type of sequence to retrieve: cds (coding sequences), pep (peptides), or both (default)."
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
invalid_out(){
    echo "Error 2: Invalid output (-o) option: $1"
    echo "Error 2: Output (-o) must be one of: both, pep, cds"
    echo "====================================="
    exit 1
}

############################################################

echo
echo "====================================="
echo "Get Ensembl sequences with rsync for all species in a file."
release="98"
input="/dev/stdin"
output="both"
# Defaults

while getopts ":r:f:o:h" arg; do
    case $arg in
        h) display_usage;;
        r) release=$OPTARG;;
        f) input=$OPTARG;;
        o) output=$OPTARG;;
        \?) invalid_opt $OPTARG ;;
    esac
done
# Parse input options

if [ $output != "both" ] && [ $output != "pep" ] && [ $output != "cds" ]; then
    invalid_out $output
fi
# Check output option

release="release-"$release
echo "Reading species names from: " $input
echo "Setting Ensembl release to: " $release
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
        rsync -av rsync://ftp.ensembl.org/ensembl/pub/$release/fasta/$spec/cds/*.cds.all.fa.gz $spec/.
        rsync -av rsync://ftp.ensembl.org/ensembl/pub/$release/fasta/$spec/cds/CHECKSUMS $spec/CHECKSUMS-CDS
    fi

    if [ $output == "cds" ] || [ $output == "both" ]; then
        rsync -av rsync://ftp.ensembl.org/ensembl/pub/$release/fasta/$spec/pep/*.pep.all.fa.gz $spec/.
        rsync -av rsync://ftp.ensembl.org/ensembl/pub/$release/fasta/$spec/pep/CHECKSUMS $spec/CHECKSUMS-PEP
    fi
    #rsync -av rsync://ftp.ensembl.org/ensembl/pub/latest/fasta/$spec/cds/ $spec/
    #rsync -av rsync://ftp.ensembl.org/ensembl/pub/latest/fasta/$spec/pep/ $spec/
done < "$input"
# Run rsync

echo
echo "Done!"
echo "====================================="
