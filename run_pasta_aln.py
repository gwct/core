#!/usr/bin/python
#############################################################################
#Makes pasta alignments from an input directory.
#
#Dependencies: pasta, core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
import core

############################################
#Function Definitions
############################################

def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Makes pasta alignments from an input directory. Dependencies: core, pasta");

	parser.add_argument("-i", dest="input_dir", help="A directory containing multiple multi-FASTA files to be aligned.");
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. -v 1: print all pasta output, -v 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output_dir", help="A directory to which the aligned sequences will be written.");

	args = parser.parse_args();

	if args.input_dir == None or args.output_dir == None:
		parser.print_help();
		sys.exit();

	if args.verbosity not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 1: -v must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	return args.input_dir, args.verbosity, args.output_dir;

############################################
#Main Block
############################################

indir, v, outdir = IO_fileParse();

if outdir[len(outdir)-1] != "/":
	outdir = outdir + "/";

print "=======================================================================";
print "Aligning all files in:\t\t" + indir;
print "Writing alignments to:\t\t" + outdir;
print "-------------------------------------";

filelist = os.listdir(indir);

if not os.path.exists(outdir):
	os.system("mkdir " + outdir);

i = 0;
numbars = 0;
donepercent = [];

for each in filelist:
	if v == 0:
		numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;
	if each.find(".fa") == -1:
		continue;

	infilename = indir + each;
	jobname = each[:each.index(".fa")];
	pastadir = outdir + jobname + "_pasta/";
	outfilename = outdir + each[:each.index(".fa")] + "_aln.fa";

	if pastadir not in os.listdir(outdir):
		mkcmd = "mkdir " + pastadir;
		os.system(mkcmd);
	
	if v == 1:
		cmd = "run_pasta.py -i " + infilename + " -d protein -j " + jobname + " -o " + pastadir;
	elif v == 0:
		cmd = "run_pasta.py -i " + infilename + " -d protein -j " + jobname + " -o " + pastadir + " > /dev/null";
	os.system(cmd);

	pastafiles = os.listdir(pastadir);
	for pfile in pastafiles:
		if pfile.find("marker") != -1:
			cpcmd = "cp " + pastadir + pfile + " " + outfilename;
			os.system(cpcmd);
			break;

if v == 0:
	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print "=======================================================================";





