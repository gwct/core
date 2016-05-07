#!/usr/bin/python
#############################################################################
# Gene family annotation lookup with GO terms.
#
# Gregg Thomas
# Feb. 2016
#############################################################################

import sys, os, platform, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Gene family annotation lookup with GO terms.");

	parser.add_argument("-i", dest="input_file", help="A tab delimited file containing on each line the gene family ID in the first column and all GO terms associated with that family.");
	parser.add_argument("-f", dest="fam_list", help="A comma delimited list of gene families for which you wish to look up GO terms. Enter 'all' to retrieve all families in input file.");
	parser.add_argument("-q", dest="query_list", help="A comma delimted list of words or terms to search for in the families retrieved with -f.", default="");
	parser.add_argument("-o", dest="output_format", help="1 - nested tabulated format. 2 - rowed tabulated format. 3 - Just print the number of families with a certain GO term (used in conjunction with -q). Default: 1", type=int, default=1);

	args = parser.parse_args();

	if errorflag == 0:
		if args.input_file == None or args.fam_list == None:
			core.errorout(1, "All options must be specified.");
			optParse(1);

		if args.output_format not in [1,2,3]:
			core.errorout(2, "-o can take values of only 1, 2, or 3.");
			optParse(1);

		if args.output_format == 3 and args.query_list == None:
			core.errorout(3, "Setting -o 3 requires you to also be searching for terms with -q.");
			optParse(1);

		return args.input_file, args.fam_list.split(","), args.query_list.split(","), args.output_format;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

######################

def printFams(fdict, of):

	#fams = sorted([int(k) for k in fdict.keys()])
	fams = fdict;

	if of == 1:
		for fam in fams:
			fam = str(fam);
			print "Family ID: " + fam + " (" + fdict[fam][0] + " genes)";
			for gt in fdict[fam]:
				if not isinstance(gt, list):
					continue;
				for g in range(len(gt)):
					if g == 1:
						continue;
					if g == 0:
						print "\t" + gt[g] + " (" + gt[g+1] + " genes)";
					else:
						print "\t\t" + gt[g];

	elif of == 2:
		for fam in fams:
			fam = str(fam);
			i = 0;
			for gt in fdict[fam]:
				if i == 0:
					i = i + 1;
					continue;
				outline = fam + "\t" + str(fdict[fam][0]) + "\t";
				#elif i > 1:
				#	outline = "\t\t";
				for g in gt:
					if any(x for x in ["name: ","namespace: ", "def: "] if x in g):
						g = g[g.index(": ")+2:]
					outline = outline + g + "\t";
				print outline;
				i = i + 1;


	#print "# ---------";
	#print "# " + str(len(fdict)) + " total families."

######################

def querySearch(fdict, qlist, ofp):
	pdict = {};
	for q in qlist:
		for fam in fdict:
			for gt in fdict[fam]:
				for g in gt:
					if any(v in g for v in [q, q.title(), q.upper(), q.lower()]):
						pdict[fam] = fdict[fam];

	if ofp in [1,2]:
		printFams(pdict, ofp);
	elif ofp == 3:
		return len(pdict);

############################################
#Main Block
############################################

infilename, famlist, querylist, outformat = optParse(0);
# Getting the input parameters.

if platform.system() == 'Windows':
	gofilename = "C:\\bin\go.obo"
else:
	gofilename = "/Users/Gregg/bin/go.obo";

# print "# =======================================================================";
# print "#\t\tGene Family Annotation Lookup";
# print "#\t\t\t" + core.getDateTime()
# print "# Input gene family annotation file:\t" + infilename;
# print "# GO term database location:\t\t" + gofilename;
# print "# ---------";
# print "# Looking up GO terms for the following", len(famlist), "families:";
# print "# " + ",".join(famlist);
# print "# ---------";
# if querylist != [""]:
# 	print "# Looking for the following phrases in above families:"
# 	for each in querylist:
# 		print "# " + each;
# 	print "# ---------";
# print "# Looking up annotations...\n"

if famlist[0].lower() == "all":
	famlist = [f.split("\t")[0] for f in open(infilename)];

# if len(famlist) > 100:
# 	if len(famlist) > 1000:
# 		print "# WARNING: Looking up more than 1000 families... this will take a very long time.\n"
# 	else:
# 		print "# WARNING: Looking up more than 100 families... this may take a while.\n"

famdict = {};
total_fams = 0;

for fam in famlist:
	famdict[fam] = [];
	famline = "";
	for line in open(infilename):
		line = line.strip().split();
		if line[0][:line[0].index("-")] == fam:
			famline = line;
			famcount = line[0][line[0].index("-")+1:];
			famdict[fam].append(famcount);
			#print line[0], famcount;
			break;


	if famline == "":
		continue;

	famline.pop(0);
	for goterm in famline:
		#print goterm;
		golist = [];
		gocount = goterm[goterm.index("-")+1:];
		goterm = goterm[:goterm.index("-")];
		golist.append(goterm);
		golist.append(gocount);
		gofile = open(gofilename, "r");	
		for line in gofile:
			if line.strip() == "id: " + goterm:
				line = gofile.next();
				while "[Term]" not in line:
					if any(l in line for l in ["name:","namespace:","def:"]):				
						golist.append(line.strip());
					line = gofile.next();
				famdict[fam].append(golist);
				break;
		gofile.close();

	if len(famdict) > 100:
		if outformat in [1,2]:
			if querylist != [""]:
				querySearch(famdict, querylist, outformat);
			else:
				printFams(famdict, outformat);
		elif outformat == 3:
			if querylist != [""]:
				total_fams = total_fams + querySearch(famdict, querylist, outformat);
			else:
				total_fams = total_fams + 100;
		famdict = {};
		#print "# Continuing lookup...\n"

if outformat in [1,2]:
	if querylist != [""]:
		querySearch(famdict, querylist, outformat);
	else:
		printFams(famdict, outformat);
elif outformat == 3:
	if querylist != [""]:
		total_fams = total_fams + querySearch(famdict, querylist, outformat);
	else:
		total_fams = total_fams + len(famdict);

# 	print "# ---------";
# 	print "# " + str(total_fams) + " total families found."

# print "\n# Done!";
# print "# =======================================================================";