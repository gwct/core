#!/usr/bin/python
#############################################################################
#Makes muscle alignments from an input directory.
#
#Dependencies: muscle, core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Makes muscle alignments from an input directory. Dependencies: core, muscle");

	parser.add_argument("-i", dest="input_dir", help="A directory containing multiple multi-FASTA files to be aligned.");
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. -v 1: print all muscle output, -v 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output_dir", help="A directory to which the aligned sequences will be written.");

	args = parser.parse_args();

	if errorflag == 0:	
		if args.input_dir == None or args.output_dir == None:
			parser.print_help();
			sys.exit();

		if args.verbosity not in [0,1]:
			core.errorOut(1, "-v must take values of either 0 or 1");
			optParse(1);

		return args.input_dir, args.verbosity, args.output_dir;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

indir, v, outdir = optParse(0);

if not os.path.isdir(indir):
	errorOut(0, "-i must be a valid directory path");
	sys.exit();
indir = os.path.abspath(indir);
if indir[-1] != "/":
	indir = indir + "/";

if outdir[len(outdir)-1] != "/":
	outdir = outdir + "/";

print("=======================================================================");
print("\t\t\tAligning with MUSCLE");
print("\t\t\t" + core.getDateTime());
print("Aligning all files in:\t\t" + indir);
print("Writing alignments to:\t\t" + outdir);
print("-------------------------------------");

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
	outfilename = outdir + each[:each.index(".fa")] + "_aln.fa";

	if v == 1:
		cmd = "muscle -in " + infilename + " -out " + outfilename;
	elif v == 0:
		cmd = "muscle -in " + infilename + " -out " + outfilename + " -quiet";
	os.system(cmd);


if v == 0:
	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
print("\n" + core.getTime() + " Done!");
print("=======================================================================");





