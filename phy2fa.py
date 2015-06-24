#!/usr/bin/python
#############################################################################
#Converts Phylip formatted sequences to FASTA format.
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

def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Converts FASTA formatted sequences to Phylip format. Dependencies: core");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many PHYLIP files or a single PHYLIP file.");
	parser.add_argument("-o", dest="output", help="Output. Either the directory where files will be written or simply an output file name.");

	args = parser.parse_args();

	if args.input == None or args.output == None:
		parser.print_help();
		sys.exit();

	return args.input, args.output;

#########

def pf_Converter(i,o):
	inseqs = core.phylipGetDict2(i);
	seqs = inseqs[0];

	outfile = open(o, "w");
	outfile.write("");
	outfile.close();

	for title in seqs:
		ftitle = ">" + title;
		core.writeSeq(o,seqs[title],ftitle);


############################################
#Main Block
############################################

ins, outs = IO_fileParse();

if os.path.isfile(ins):
	print "Input type is:\tFile"
	print "Converting to FASTA format..."

	pf_Converter(ins,outs);

	print "Done!";

else:
	print "Input type is:\tDirectory"
	print "Converting all Phylip files to FASTA format..."

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
		if each.find(".phy") == -1:
			continue;

		infilename = ins + each;
		outfilename = outs + each[:each.index(".phy")] + ".fa";

		pf_Converter(infilename,outfilename);

	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
	print "\nDone!";


