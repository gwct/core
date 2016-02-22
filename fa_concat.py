#!/usr/bin/python
#############################################################################
#This script concatenates FASTA formatted sequences in multiple files to a single
#FASTA file.
#
#Usage: python count_aln.py [input directory] [output file name]
#
#Depenencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

if len(sys.argv) != 3 or sys.argv[1] == '-h':
	print "Usage:\t$ fa_conact.py [input directory] [output file name]";
	sys.exit();

ins = sys.argv[1];
outs = sys.argv[2];

print "=======================================================================";
print "Concatenating alignments in:\t\t" + ins;
print "Writing concatenated alignment to:\t" + outs;
print "-------------------------------------";

filelist = os.listdir(ins);

if outs.find("/") != -1:
	pos = outs.rfind('/');
	outdir = outs[:pos] + '/';
	if not os.path.exists(outdir):
		print "+Creating output directory.";
		os.system("mkdir " + outdir);
	partfile = outdir + "partitions.txt";
else:
	partfile = "partitions.txt";

pfile = open(partfile, "w");
pfile.write("");

concats = {};
numpos = 0;

i = 0;
numbars = 0;
donepercent = [];

for each in filelist:
	numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;

	if each.find(".fa") == -1:
		continue;

	infilename = ins + each;

	inseqs = core.fastaGetDict(infilename);
	seqlen = len(inseqs[inseqs.keys()[0]]);

	pfile.write(each[:each.index("_")] + " = " + str(numpos+1) + "-" + str(numpos+seqlen) + "\n");
	numpos = numpos + seqlen;

	for title in inseqs:
		if " " in title:
			newtitle = title[:title.index(" ")];
		else:
			newtitle = title;
		if newtitle not in concats:
			concats[newtitle] = inseqs[title];
		else:
			concats[newtitle] = concats[newtitle] + inseqs[title];

pfile.close();
outfile = open(outs, "w");
for spec in concats:
	outfile.write(spec);
	outfile.write("\n");
	outfile.write(concats[spec]);
	outfile.write("\n");
outfile.close();

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print "=======================================================================";
print numpos;
