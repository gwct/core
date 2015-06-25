#!/usr/bin/python
#############################################################################
#Converts FASTA formatted sequences to Phylip format.
#
#Dependencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################

def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Converts FASTA formatted sequences to Phylip format. Dependencies: core");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-o", dest="output", help="Output. Either the directory where files will be written or simply an output file name.");

	args = parser.parse_args();

	if args.input == None or args.output == None:
		parser.print_help();
		sys.exit();

	return args.input, args.output;

#########

def fp_Converter(i,o):
	inseqs = core.fastaGetDict(i);

	outfile = open(o, "w");
	outline = " " + str(len(inseqs)) + " " + str(len(inseqs[inseqs.keys()[0]])) + "\n";
	outfile.write(outline);

	interval = len(max(inseqs.keys(), key=len)) + 3;

	for title in inseqs:
		newtitle = title[1:];
		spaces = " " * (interval - len(newtitle));
		seq = inseqs[title];

		outline = newtitle + spaces + seq + "\n";
		outfile.write(outline);

	outfile.close();

############################################
#Main Block
############################################

ins, outs = IO_fileParse();

if os.path.isfile(ins):
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "Input type is:\tFile"
	print "Converting to Phylip format..."

	fp_Converter(ins,outs);

	print core.getTime() + " Done!";
	print "=======================================================================";

else:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "Input type is:\tDirectory"
	print "Converting all FASTA files to Phylip format..."

	if not ins.endswith("/"):
		ins = ins + "/";
	if not outs.endswith("/"):
		outs = outs + "/";

	if not os.path.exists(outs):
		os.system("mkdir " + outs);

	filelist = os.listdir(ins);

	numfiles = len(filelist);
	numbars = 0;
	donepercent = [];
	k = 0;

	for each in filelist:
		numbars, donepercent = core.loadingBar(k, numfiles, donepercent, numbars);
		k = k + 1;
		if each.find(".fa") == -1:
			continue;

		infilename = ins + each;
		outfilename = outs + each[:each.index(".fa")] + ".phy";

		fp_Converter(infilename,outfilename);

	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	print "=======================================================================";

