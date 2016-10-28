#!/usr/bin/python
#############################################################################
#Script to count the total number of positions in a fasta file or a directory full of fasta files.
#
#Usage: python count_pos.py [input file or directory] [1,0]
#
#The script first checks if the input is a file or directory. If it is a file it will just count the
#positions in that file and display. If it is a directory it will count the number of positions in
#all files and print the sum. If the second parameter is set to 1, it will also print the number of
#positions in each file separately.
#
#Dependencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

if len(sys.argv) not in [1,2,3]:
	print "Usage:\t$ count_pos.py [input directory or filename] [1,0 to display individual file counts or not]";
	sys.exit();

ins = sys.argv[1];
disp_file = 0;
if len(sys.argv) > 2:
	disp_file = sys.argv[2];
if disp_file not in ["0","1"]:
	print "Not printing file counts.";
	disp_file = 0;

disp_file = int(disp_file);

print "=======================================================================";
print "\t\t\t" + core.getDateTime();
print "Counting the total number of positions (AAs or NTs) in:\t" + ins;

if os.path.isfile(ins):
	if disp_file == 1:
		print "----------";
		print "Sequence\tLength";
	inseqs = core.fastaGetDict(ins);
	tot_pos = 0;
	for seq in inseqs:
		if disp_file == 1:
			print seq + "\t" + str(len(inseqs[seq]));
		tot_pos = tot_pos + len(inseqs[seq]);
	print "----------";
	print "Total sequences:\t" + str(len(inseqs));
	print "Total positions:\t" + str(tot_pos);
	print "=======================================================================";

else:
	if not ins.endswith("/"):
		ins = ins + "/";
	filelist = os.listdir(ins);

	tot_pos = 0;

	numlines = len(filelist);
	numbars = 0;
	donepercent = [];
	i = 0;

	for each in filelist:

		if disp_file == 0:
			numbars, donepercent = core.loadingBar(i, numlines, donepercent, numbars);
		elif disp_file == 1:
			print "----------";
			print each;
			print "Sequence\tLength";
		i = i + 1;

		if each.find(".fa") == -1:
			continue;

		specpos = 0;

		infilename = ins + each;

		inseqs = core.fastaGetDict(infilename);

		for seq in inseqs:
			tot_pos = tot_pos + len(inseqs[seq]);
			if disp_file == 1:
				specpos = specpos + len(inseqs[seq]);
				print seq + "\t" + str(len(inseqs[seq]));

		if disp_file == 1:
			print "Total\t" + str(specpos);

	if disp_file == 0:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	elif disp_file == 1:
		print "----------";
	print "\n" + core.getTime() + " Done!";
	print "-----";
	print "Total residues:\t", tot_pos;
	print "=======================================================================";

