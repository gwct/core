#!/usr/bin/python
############################################
#CAFE helper script to give some very simple
#stats about the input file.
#
#Sample usage: python cafe_get_stats.py [input_file_name]
############################################

import sys

infilename = sys.argv[1];

i = 0;
numgenes = 0;

for line in open(infilename):
	tmpline = line.replace("\n","").split("\t");
	if i == 0:
		numspecs = len(tmpline)-2;

		i = i + 1;
		continue;
	
	for x in xrange(len(tmpline)):
		if x <= 1:
			continue;
		else:
			numgenes = numgenes + int(tmpline[x]);

	i = i + 1;

print "====================================";
print infilename;
print "Number of species:\t", numspecs;
print "Number of families:\t", i-1;
print "Number of genes:\t", numgenes;
print "====================================";
