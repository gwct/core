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
#############################################################################
import sys, os
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

ins = sys.argv[1];
print "=======================================================================";
print "\t\t\t" + core.getDateTime();
print "Counting the total number of positions (AAs or NTs) in:\t" + ins;

if os.path.isfile(ins):
	inseqs = core.fastaGetDict(ins);
	tot_pos = 0;
	for seq in inseqs:
		tot_pos = tot_pos + len(inseqs[seq]);
	print ins + "\t" + str(tot_pos);

else:
	filelist = os.listdir(ins);
	disp_file = 0;
	if len(sys.argv) > 2:
		disp_file = sys.argv[2];
	if disp_file not in ["0","1"]:
		print "Not printing file counts.";
		disp_file = 0;

	disp_file = int(disp_file);
	tot_pos = 0;

	numlines = len(filelist);
	numbars = 0;
	donepercent = [];
	i = 0;

	for each in filelist:

		if disp_file == 0:
			numbars, donepercent = core.loadingBar(i, numlines, donepercent, numbars);
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

		if disp_file == 1:
			print each + "\t" + str(specpos);

	if disp_file == 0:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	print "-----";
	print "Total residues:\t", tot_pos;
	print "=======================================================================";

