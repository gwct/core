#!/usr/bin/python
#############################################################################
#Gets the output from a run_raxml.py (gene trees and bootstrap replicate 
#files) run ready for input into ASTRAL.
#
#Dependencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################
import sys, os
import core

if len(sys.argv) > 1:
	indir = os.path.abspath(sys.argv[1]) + "/";
else:
	indir = os.getcwd() + "/";
raxtreedir = indir + "raxml_best/";
raxoutdir = indir + "raxml_out/";
outfilename = indir + "astral_gt.txt";
bsfilename = indir + "astral_bs.txt";

print "=======================================================================";
print "\tGetting bootstrap replicates for each gene\n\ttree from a run_raxml.py run for input to ASTRAL";
print "\t\t" + core.getDateTime()
print "INPUT    | Input directory:\t\t\t" + indir;
print "INFO     | All other directories and files will be within the input directory.";
print "INFO     | RAxML best tree directory:\t\traxml_best/";
print "INFO     | RAxML output directory:\t\traxml_out/";
print "OUTPUT   | Writing genetrees to:\t\tastral_gt.txt";
print "OUTPUT   | Writing locations of BS files to:\tastral_bs.txt";
print "-------------------------------------";
print core.getTime() + " | Working...";

ofile = open(outfilename, "w");
ofile.write("");
bfile = open(bsfilename, "w");
bfile.write("");

treelist = os.listdir(raxtreedir);

r = 0;
d = 10;
i = 0;

for each in treelist:
	i = i + 1;
	core.loadingRotator(i,r,d);
	if each.find("RAxML_bestTree.") == -1:
		continue;
	gid = each[each.index(".")+1:];

	gtree = open(raxtreedir + each,"r").read();
	gbsfile = raxoutdir + gid + "_raxout/RAxML_bootstrap." + gid;

	ofile.write(gtree)
	bfile.write(gbsfile + "\n");

ofile.close();
bfile.close();
sys.stderr.write('\b');
print core.getTime() + " | Done!"
print "=======================================================================";
