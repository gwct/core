#!/usr/bin/python
#############################################################################
#Script to gather counts and info about alignments.
#
#Usage: python count_aln.py [input file or directory] [1,0]
#
#The script first checks if the input is a file or directory. If it is a directory it will gather info in
#all files and print the sums. If the second parameter is set to 1, it will also print the number of
#positions in each file separately.
#
#Depenencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################
def alignCounter(ifile):
	#print ifile;
	inseqs = core.fastaGetDict(ifile);
	num_seqs = len(inseqs);
	tot_pos = 0;
	seq_len = 0;

	inv_sites = 0;
	var_sites = 0;
	gaps = 0;
	gap_sites = 0;
	gap_dist = {};
	for x in range(num_seqs):
		gap_dist[x+1] = 0;

	j = 0;

	for seq in inseqs:
		if j == 0:
			j = j + 1;
			seq_len = len(inseqs[seq]);
		tot_pos = tot_pos + len(inseqs[seq]);

	for col in range(seq_len):
		site = [];
		for seq in inseqs:
			site.append(inseqs[seq][col]);

		gaps = gaps + site.count("-");

		if site.count("-") != len(site):
			for base in site:
				if base != "-":
					if site.count(base) == (len(site) - site.count("-")):
						inv_sites = inv_sites + 1;
					else:
						var_sites = var_sites + 1;
					break;

		if "-" in site:
			gap_sites = gap_sites + 1;
			gap_dist[site.count("-")] = gap_dist[site.count("-")] + 1;

	if disp_file == 1:
		print(ifile + "\t" + str(num_seqs) + "\t" + str(tot_pos) + "\t" + str(seq_len) + "\t" + str(inv_sites) + "\t" + str(var_sites) + "\t" + str(gaps) + "\t" + str(gap_sites) + "\t" + str(gap_dist));
	return num_seqs,tot_pos,seq_len,inv_sites,var_sites,gaps,gap_sites,gap_dist;

############################################
#Main Block
############################################
if len(sys.argv) not in [1,2,3]:
	print("Usage:\t$ count_aln.py [input directory or filename] [1,0 to display individual file counts or not]");
	sys.exit();

ins = sys.argv[1];

if os.path.isfile(ins):
	print("=======================================================================");
	print("\t\t\t" + core.getDateTime());
	print("Gathering alignment stats from file:\t" + ins);
	print("-----");
	disp_file = 1;
	print("Filename\t# Seqs\tTotal # Positions\tAlignment Length\t# Invariant Sites\t# Variant Sites\t# Gaps\t# Sites with Gaps\tGap Distribution");
	numseqs, totpos, seqlen, invsites, varsites, gap, gapsites, gapdist = alignCounter(ins);
	print(core.getTime() + " Done!");
	print("=======================================================================");

else:
	print("=======================================================================");
	print("\t\t\t" + core.getDateTime());
	print("Gathering alignment stats from all .fa files in:\t" + ins);
	filelist = os.listdir(ins);
	disp_file = 0;
	if len(sys.argv) > 2:
		disp_file = sys.argv[2];
	if disp_file not in ["0","1"]:
		print("Not printing file counts.");
		disp_file = 0;

	disp_file = int(disp_file);

	numlines = len(filelist);
	numbars = 0;
	donepercent = [];
	i = 0;

	totseqs = 0;
	allpos = 0;
	alllen = 0;
	totinv = 0;
	totvar = 0;
	totgap = 0;
	totgapsite = 0;

	for each in filelist:
		if disp_file == 0:
			numbars, donepercent = core.loadingBar(i, numlines, donepercent, numbars);
		i = i + 1;

		if each.find(".fa") == -1:
			continue;

		infilename = ins + each;

		numseqs, totpos, seqlen, invsites, varsites, gap, gapsites, gapdist = alignCounter(infilename);

		totseqs = totseqs + numseqs;
		allpos = allpos + totpos;
		alllen = alllen + seqlen;
		totinv = totinv + invsites;
		totvar = totvar + varsites;
		totgap = totgap + gap;
		totgapsite = totgapsite + gapsites;
		
	if disp_file == 0:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print("\n" + core.getTime() + " Done!");
	#print "# Seqs\tTotal # Positions\tTotal # Columns\t# Invariant Sites\t# Variant Sites\t# Gaps\t# Sites with Gaps";
	#print str(totseqs) + "\t" + str(allpos) + "\t" + str(alllen) + "\t" + str(totinv) + "\t" + str(totvar) + "\t" + str(totgap) + "\t" + str(totgapsite);
	print("-----");
	print("Total # Seqs\t" + str(totseqs));
	print("Total # Positions\t" + str(allpos));
	print("Total # Columns\t" + str(alllen));
	print("# Invariant Sites\t" + str(totinv));
	print("# Variant Sites\t" + str(totvar));
	print("Total # Gaps\t" + str(totgap));
	print("# Sites with Gaps\t" + str(totgapsite));
	print("=======================================================================");



