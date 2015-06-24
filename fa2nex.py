#!/usr/bin/python
#############################################################################
#Converts FASTA formatted sequences to Nexus format.
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

	parser = argparse.ArgumentParser(description="Converts FASTA formatted sequences to Nexus format. Dependencies: core");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-o", dest="output", help="Output. A single file to which the (concatenated) sequences will be written in NEXUS format.");

	args = parser.parse_args();

	if args.input == None or args.output == None:
		parser.print_help();
		sys.exit();

	return args.input, args.output;


############################################
#Main Block
############################################

ins, outs = IO_fileParse();

outfile = open(outs, "w");
outfile.write("#nexus\n\nbegin data;\n");
outfile.close();

seqs = {};

if os.path.isfile(ins):
	print "Input type is:\tFile"
	print "Converting to Nexus format..."

	inseqs = core.fastaGetDict(ins);

	for title in inseqs:
		if title in seqs:
			seqs[title] = seqs[title] + inseqs[title].replace(" ","");
		else:
			seqs[title] = inseqs[title].replace(" ","");

else:
	print "Input type is:\tDirectory"
	print "Reading all FASTA sequences..."

	if not ins.endswith("/"):
		ins = ins + "/";

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
		inseqs = core.fastaGetDict(infilename);

		for title in inseqs:
			if title in seqs:
				seqs[title] = seqs[title] + inseqs[title].replace(" ","");
			else:
				seqs[title] = inseqs[title].replace(" ","");
		

	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);

print "Writing concatenated sequences to nexus file...";

ntaxa = len(seqs);
nsites = len(seqs[seqs.keys()[0]]);

outfile = open(outs, "a");
outline = "\tdimensions ntax=" + str(ntaxa) + " nchar=" + str(nsites) + ";\n";
outfile.write(outline);
outline = "\tformat datatype=protein gap=-;\n";
outfile.write(outline);
outline = "\tmatrix\n";
outfile.write(outline);

interval = len(max(seqs.keys(), key=len)) + 3;

for title in seqs:
	newtitle = title[1:];
	spaces = " " * (interval - len(newtitle));
	seq = seqs[title];

	outline = newtitle + spaces + seq + "\n";
	outfile.write(outline);

outline = "\t;\n";
outfile.write(outline);
outline = "end;";
outfile.write(outline);

outfile.close();
print "\nDone!";



