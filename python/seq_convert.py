#!/usr/bin/python
#############################################################################
#Converts sequence formats from one to another. Converts between any of
#FASTA, Phylip, and Nexus formats. Please note, this script assumes the
#file extensions of .fa, .ph, and .nex, respectively, for those formats.
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

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Converts sequence formats from one to another. Converts between any of FASTA, Phylip, and Nexus formats. Please note, this script assumes the file extensions of .fa, .ph, and .nex, respectively, for those formats. Dependencies: core");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many sequence files or a single sequence file.");
	parser.add_argument("-f", dest="input_type", help="The format of the input sequences.");
	parser.add_argument("-o", dest="output", help="Output. Either the directory where files will be written or simply an output file name.");
	parser.add_argument("-t", dest="output_type", help="The desired output format of the sequences.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.input == None or args.output == None:
			parser.print_help();
			sys.exit();

		intype = args.input_type.lower();
		outtype = args.output_type.lower();

		for t in [intype, outtype]:
			if t not in ["fasta", "phylip", "nexus", "fa", "phy", "ph", "nex", "f", "p", "n"]:
				core.errorOut(1, "-f and -t must take values of fasta, nexus, or phylip");
				optParse(1);

		return args.input, intype[:1], args.output, outtype[:1];

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

ins, fr, outs, to = optParse(0);

if os.path.isfile(ins):
	fileflag = 1;
	filelist = [ins];
else:
	fileflag = 0;
	if not os.path.isdir(ins):
		errorOut(0, "-i must be a valid directory path");
		sys.exit();
	ins = os.path.abspath(ins);
	if ins[-1] != "/":
		ins = ins + "/";
	outs = os.path.abspath(outs);
	if not outs.endswith("/"):
		outs = outs + "/";
	filelist = os.listdir(ins);

print("==============================================================================================");
print("\t\t\tSequence format conversion");
print("\t\t\t" + core.getDateTime());
if fileflag == 1:
	print("INPUT    | Converting file: " + ins);
else:
	print("INPUT    | Converting all files from directory: " + ins);
print("INFO     | Input format:  " + fr);
print("INFO     | Output format: " + to);
if fileflag == 1:
	print("OUTPUT   | Writing output to file: " + outs);
else:
	print("OUTPUT   | Writing output files to directory:   " + outs);
print("-------------------------------------");

if fileflag == 0:
	print(core.getTime() + " | Creating output directory...");
	if not os.path.exists(outs):
		os.system("mkdir " + outs);

	numfiles = len(filelist);
	numbars = 0;
	donepercent = [];
	i = 0;

if fr == "f":
	init = ".fa";
elif fr == "p":
	init = ".ph";
elif fr == "n":
	init = ".nex";

if to == "f":
	suffix = ".fa";
elif to == "p":
	suffix = ".ph";
elif to == "n":
	suffix = ".nex";

firstbar = True

for each in filelist:
	if fileflag == 0:
		numbars, donepercent, firstbar = core.loadingBar(i, numfiles, donepercent, numbars, firstbar=firstbar, disperc=True);
		i = i + 1;

	if fr == "f" and each.find(".fa") == -1:
		continue;
	if fr == "p" and each.find(".ph") == -1:
		continue;
	if fr == "n" and each.find(".nex") == -1:
		continue;

	if fileflag == 1:
		infilename = each;
		outfilename = outs;
	else:
		infilename = ins + each;
		outfilename = outs + each[:each.index(init)] + suffix;

	outfile = open(outfilename, "w");
	outfile.write("");
	outfile.close();
	
	if fr == "f":
		inseqs = core.fastaGetDict(infilename);
		newseqs = {};
		for title in inseqs:
			newtitle = title[1:];
			newseqs[newtitle] = inseqs[title];
		inseqs = newseqs;
	if fr == "p":
		inseqs = core.phylipGetDict(infilename)[0];
	if fr == "n":
		inseqs = core.nexusGetDict(infilename);

	newseqs = {};
	for title in inseqs:
		newtitle = title.replace(" ",":");
		newseqs[newtitle] = inseqs[title];
	inseqs = newseqs;

	if to == "f":
		newseqs = {};
		for title in inseqs:
			newtitle = ">" + title;
			newseqs[newtitle] = inseqs[title];
		inseqs = newseqs;
		core.writeFasta(inseqs,outfilename);
	if to == "p":
		core.writePhylip(inseqs,outfilename);
	if to == "n":
		core.writeNexus(inseqs,outfilename);

if fileflag == 0:
	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
	print("\n");

print(core.getTime() + " | Done!");
print("==============================================================================================");
