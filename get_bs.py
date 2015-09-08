#!/usr/bin/python
#############################################################################
#Gets the output from a run_raxml.py (gene trees and bootstrap replicate 
#files) run ready for input into ASTRAL.
#
#Dependencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################
import sys, os, argparse
import core

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Compiles the info from a run_raxml run into files used by ASTRAL (tree file and bootstrap locations file). Input directory must be a run_raxml output directory.");

	parser.add_argument("-i", dest="input", help="A run_raxml output directory. Default: current directory.", default=os.getcwd()+"/");
	parser.add_argument("-b", dest="bs_relabel", help="Set to 0 to get original RAxML bootstrap files. Set to 1 to get relabeled files. Default: 0", type=int, default=0);

	args = parser.parse_args();

	if errorflag == 0:

		if args.bs_relabel not in [0,1]:
			core.errorOut(1, "-r must take values of either 0 or 1");
			optParse(1);
	
		return args.input, args.bs_relabel;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();


############################################
#Main Block
############################################
indir, bs_opt = optParse(0);

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
if bs_opt == 1:
	print "INFO     | Retrieving bootstraps from relabeled files.";
elif bs_opt == 0:
	print "INFO     | Retrieving bootstraps from original RAxML files.";
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
	if bs_opt == 0:
		gbsfile = raxoutdir + gid + "_raxout/RAxML_bootstrap." + gid;
	elif bs_opt == 1:
		gbsfile = raxoutdir + gid + "_raxout/bootstrap_relabel.txt";

	ofile.write(gtree)
	bfile.write(gbsfile + "\n");

ofile.close();
bfile.close();
sys.stderr.write('\b');
print core.getTime() + " | Done!"
print "=======================================================================";
