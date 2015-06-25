#!/usr/bin/python
#############################################################################
#Converts FASTA formatted sequences to Nexus format.
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

	parser = argparse.ArgumentParser(description="Converts FASTA formatted sequences to Nexus format. Dependencies: core");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-o", dest="output", help="Output. An output directory or filename to which the sequences will be written in NEXUS format.");

	args = parser.parse_args();

	if args.input == None or args.output == None:
		parser.print_help();
		sys.exit();

	return args.input, args.output;

#####
def writeNexus(seqdict, ofname):
	outfile = open(ofname, "w");
	outfile.write("#NEXUS\n\nBegin data;\n");

	ntaxa = len(seqdict);
	nsites = len(seqdict[seqdict.keys()[0]]);

	outline = "\tDimensions ntax=" + str(ntaxa) + " nchar=" + str(nsites) + ";\n";
	outfile.write(outline);
	outline = "\tFormat datatype=protein gap=-;\n";
	outfile.write(outline);
	outline = "\tMatrix\n";
	outfile.write(outline);

	interval = len(max(seqdict.keys(), key=len)) + 3;

	for title in seqdict:
		newtitle = title[1:];
		spaces = " " * (interval - len(newtitle));
		seq = seqdict[title];

		outline = newtitle + spaces + seq + "\n";
		outfile.write(outline);

	outline = "\t;\n";
	outfile.write(outline);
	outline = "End;";
	outfile.write(outline);

	outfile.close();

############################################
#Main Block
############################################

ins, outs = IO_fileParse();

seqs = {};

if os.path.isfile(ins):
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "Input type is:\tFile"
	print "Converting to Nexus format..."

	inseqs = core.fastaGetDict(ins);

	for title in inseqs:
	#	if title in seqs:
	#		seqs[title] = seqs[title] + inseqs[title].replace(" ","");
	#	else:
		seqs[title.replace(" ",":")] = inseqs[title].replace(" ","");
	writeNexus(seqs,  outs);

else:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();	
	print "Input type is:\tDirectory"
	print "Converting to FASTA format..."

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
		seqs = {};
		numbars, donepercent = core.loadingBar(k, numfiles, donepercent, numbars);
		k = k + 1;
		if each.find(".fa") == -1:
			continue;

		infilename = ins + each;
		outfilename = outs + each[:each.index(".fa")] + ".nex";
		inseqs = core.fastaGetDict(infilename);

		for title in inseqs:
		#	if title in seqs:
		#		seqs[title] = seqs[title] + inseqs[title].replace(" ","");
		#	else:
			seqs[title.replace(" ",":")] = inseqs[title].replace(" ","");
		writeNexus(seqs, outfilename);		


	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);

#print "Writing concatenated sequences to nexus file...";
print "\n" + core.getTime() + " Done!";
print "=======================================================================";



