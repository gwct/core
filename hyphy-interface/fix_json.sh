#!/bin/bash
#PURPOSE: Replace 'inf' with 'Infinity' and 'null' with 'NaN'; this is necessary to make HyPhy output compatible with Python

if [ "$#" -ne 1 ]; then
    echo "Error: an input file command line argument is required"
    exit 1
fi

infile=$1
echo "Input file: ${infile}"
#Comma is so that words like "information" don't get replaced
sed -i 's/inf,/Infinity,/g' ${infile}
sed -i 's/null/NaN/g' ${infile}
echo "Done!"
