#!/usr/bin/python
#############################################################################
# A script to root a set of gene trees at a given outgroup.
# 
# Dependencies: Newick Utilities
#
# Gregg Thomas, Summer 2017
#############################################################################

import sys, os

if len(sys.argv) != 4 or "-h" in sys.argv:
	print "\nUsage:\t$ python fotc.py [tree file] [outgroup] [output file name]";
	print "---->\t[tree file]: A file containing one or more unrooted, newick formatted trees, one per line.";
	print "---->\t[outgroup]: A tip label in all present in all trees at which they will be rooted.";
	print "---->\t[output filename]: The name of a file to which the rooted trees will be written.";
	sys.exit();

treefilename = sys.argv[1];
outgroup = sys.argv[2];
outfilename = sys.argv[3];

if not os.path.isfile(treefilename):
	sys.exit("\n** Error! Invalid file name for your tree file!\n");

tmpnum = 8888;
tmpfile1 = "tmp" + str(tmpnum);
while os.path.isfile(tmpfile1):
	tmpnum += 1;
	tmpfile1 = "tmp" + str(tmpnum);
tmpnum = 9839;
tmpfile2 = "tmp" + str(tmpnum);
while os.path.isfile(tmpfile2):
	tmpnum += 1;
	tmpfile2 = "tmp" + str(tmpnum);

outfile = open(outfilename, "w");
for line in open(treefilename):
	tree = line.strip();

	tfile = open(tmpfile1, "w");
	tfile.write(tree);
	tfile.close();

	cmd = "nw_reroot " + tmpfile1 + " " + outgroup + " > " + tmpfile2;
	os.system(cmd);

	pruned_tree = open(tmpfile2, "r").read();

	outfile.write(pruned_tree);
outfile.close();
os.system("rm " + tmpfile1);
os.system("rm " + tmpfile2);

