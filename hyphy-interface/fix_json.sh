#!/bin/bash
#PURPOSE: Replace 'inf' with 'Infinity'; this is necessary because HyPhy outputs infinity values in json files as 'inf' or '-inf' but Python expects 'Infinity' or '-Inifinity'

if [ "$#" -ne 1 ]; then
    echo "Error: an input file command line argument is required"
    exit 1
fi

infile=$1
echo "Input file: ${infile}"
#Comma is so that words like "information" don't get replaced
sed -i 's/inf,/Infinity,/g' ${infile}
echo "Done!"
