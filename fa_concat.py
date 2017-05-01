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

if len(sys.argv) not in [3,4] or sys.argv[1] == '-h':
	print "Usage:\t$ fa_conact.py [input directory] [output file name] [species delimiter (optional)]";
	print "For a delimiter of ' ', enter the word 'space' for that option.";
	sys.exit();

ins = sys.argv[1];
if not os.path.isdir(ins):
	print "*** Error: Input directory not found.";
outs = sys.argv[2];

print "=======================================================================";
print "Concatenating alignments in:\t\t" + ins;
print "Writing concatenated alignments to:\t" + outs;
if len(sys.argv) == 4:
	delim = sys.argv[3];
	print "Species labels should be at beginning of sequence title."
	print "Species delimiter:\t\t\t" + delim;
	print "Pre-reading all files to fill in missing data.";
	if delim.lower() == 'space':
		delim = ' ';
print "-------------------------------------";

filelist = os.listdir(ins);
if outs.find("/") != -1:
	pos = outs.rfind('/');
	outdir = outs[:pos] + '/';
	if not os.path.exists(outdir):
		print "+Creating output directory.";
		os.system("mkdir " + outdir);
	partfile = os.path.join(outdir, "partitions.txt");
else:
	partfile = "partitions.txt";

specs = [];
concats = {};
for seqfile in filelist:
	if seqfile.find(".fa") == -1:
		continue;
	infilename = os.path.join(ins, seqfile);
	inseqs = core.fastaGetDict(infilename);
	for title in inseqs:
		cur_spec = title[1:title.index(delim)];
		if cur_spec not in specs:
			specs.append(cur_spec);
		if cur_spec not in concats:
			concats[">" + cur_spec] = "";

pfile = open(partfile, "w");
pfile.write("");

numpos = 0;

i = 0;
numbars = 0;
donepercent = [];

for seqfile in filelist:
	numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;
	if seqfile.find(".fa") == -1:
		continue;
	infilename = os.path.join(ins, seqfile);
	inseqs = core.fastaGetDict(infilename);
	seqlen = len(inseqs[inseqs.keys()[0]]);

	pfile.write(seqfile[:seqfile.index("_")] + " = " + str(numpos+1) + "-" + str(numpos+seqlen) + "\n");
	numpos = numpos + seqlen;

	cur_specs = [];
	for title in inseqs:
		cur_spec = title[1:title.index(delim)];
		cur_specs.append(cur_spec);
		concats[">" + cur_spec] += inseqs[title];

	for spec in specs:
		if spec not in cur_specs:
			concats[">" + spec] += "N" * seqlen;

	# for title in inseqs:
	# 	if " " in title:
	# 		newtitle = title[:title.index(" ")];
	# 	else:
	# 		newtitle = title;
	# 	if newtitle not in concats:
	# 		concats[newtitle] = inseqs[title];
	# 	else:
	# 		concats[newtitle] = concats[newtitle] + inseqs[title];

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
