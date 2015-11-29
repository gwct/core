#!/usr/bin/python
#############################################################################
#This program reads the CAFE report file and gives the relevant results in
#a more readable format.
#
#This one has a lot of options, mostly dealing with the format you want to write the
#changing families from your species of interest. Briefly, you can specify your
#species of interest (ie "Cat") with -s, specify its ingroups and outgroups with
#-g (see below for formatting), and then write the changing families from that
#species in either a nice and readable list (-o), CAFE input format (-c, -t), or
#MCL format (-m, -u). Be sure to read through all the options below.
#
#Sample usage: python cafe_report_analysis.py -i cat_mcl_postreport.cafe -s Cat -e 1 -o cat_rapids.txt
#
#Written by: Gregg Thomas - September, 2013
#############################################################################

import sys, os, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../corelib/"))
import core

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="report_input_file", help="A CAFE report file (.cafe).");
	parser.add_argument("-s", dest="spec_of_interest", help="The species label of a particular species of interest to write families to output and estimate how many gene families have been lost in that species.", default="");
	parser.add_argument("-p", dest="output_mode", help="Option to write the output table for all species (0) or just the species of interest (1). Also, if -p 1 is specified with an MCL input/output file, only the genes from the species of interest will be extracted. Default: 0", type=int, default=0);
	parser.add_argument("-e", dest="extract_mode", help="Extract mode specifies whether only families in which the species of interest is rapidly evolving should be written to the output file (1) or if all families in which the species of interest is changing should be written (2). Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="report_output_file", help="The name of the output file. With -o, output will be formatted as in the CAFE report file.");
	parser.add_argument("-m", dest="mcl_input_file", help="The name of the original MCL dump file from which the gene family data was converted. This is to obtain the list of genes for each family in which the species of interest is changing.");
	parser.add_argument("-u", dest="mcl_output_file", help="A name for output to be written in MCL format. To be used in conjunction with -m and -e.");
	parser.add_argument("-c", dest="cafe_input_file", help="The name of the CAFE input file to obtain only families in which the species of interest is changing in CAFE format.");
	parser.add_argument("-t", dest="cafe_output_file", help="A name for output to be written in CAFE format. To be used in conjunction with -c and -t.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.report_input_file == None:
			core.errorOut(1, "A CAFE report file must be specified with -i");
			optParse(1);


		if args.output_mode not in [0,1]:
			core.errorOut(2, " -p must take values of either 0 or 1");
			optParse(1);

		if args.output_mode == 1 and args.spec_of_interest == "":
			core.errorOut(3, "-s must be specified with -p 1");
			optParse(1);

		if args.extract_mode not in [1,2]:
			core.errorOut(4, "-e must take values of either 1 or 2");
			optParse(1);

		if args.mcl_input_file != None and args.mcl_output_file == None:
			core.errorOut(5, "An MCL input file has been specified (-m) without an MCL output file (-u). Always set -u with -m.");
			optParse(1);

		if args.cafe_input_file != None and args.cafe_output_file == None:
			core.errorOut(6, "A CAFE input file has been specified (-c) without a CAFE output file (-t). Always set -t with -c.");

		return args.report_input_file, args.spec_of_interest, args.extract_mode, args.report_output_file, args.mcl_input_file, args.mcl_output_file, args.cafe_input_file, args.cafe_output_file, args.output_mode;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################

def getBranchLength(nodeid, tree):
#This function takes a node id and a tree and returns the branch length of that
#node within the tree.

	stopopen = tree.find("(", tree.index(nodeid))
	stopclose = tree.find(")", tree.index(nodeid))
	stopcomma = tree.find(",", tree.index(nodeid))

	if stopopen == -1:
		stopopen = 10000000;
	if stopclose == -1:
		stopclose = 10000000;
	if stopcomma == -1:
		stopcomma = 10000000;

	curstop = min(stopopen, stopclose, stopcomma);

	branch = tree[tree.index(":", tree.index(nodeid))+1:curstop];

	return branch;

############################################
def writeRapids(In, Out, expands, contracts, specofint, cafeinfile, em, nume, numc, etype, xtra):
#This function takes two lists of families (expands and contracts) and writes them in that format
#to an output file.

	eoutFile = open(Out, "w");
	if etype == "report":
		outline = "# Species of Interest: " + specofint + "\n# CAFE report file name: " + cafeinfile;
	elif etype == "mcl":
		outline = "# Species of Interest: " + specofint + "\n# CAFE report file name: " + cafeinfile + "\n# MCL dump file name: " + In;
	elif etype == "cafe":
		outline = "# Species of Interest: " + specofint + "\n# CAFE report file name: " + cafeinfile + "\n# CAFE input file name: " + In;
	eoutFile.write(outline);
	eoutFile.write("\n\n");

	if etype != "mcl":
		eoutFile.write("\n##################################\n");
		if em == 1:
			eoutFile.write("# Rapid Expansions in ");
			eoutFile.write(specofint)
			eoutFile.write(" (");
			eoutFile.write(str(nume));
			eoutFile.write(")");
		elif em == 2:
			eoutFile.write("# ALL Expansions in ");
			eoutFile.write(specofint)
			eoutFile.write(" (");
			eoutFile.write(str(nume));
			eoutFile.write(")");
		eoutFile.write("\n##################################\n");

	if xtra != "":
		eoutFile.write(xtra);

	for fam in expands:
		if etype == "mcl":
			fam = fam.replace("\n","").split("\t");
			for x in xrange(len(fam)):
				if x == 0:
					outline = fam[x] + " (" + str(len(fam)-1) + " genes) +\n";
					eoutFile.write(outline);
				else:
					outline = fam[x] + "\n";
					eoutFile.write(outline);
			eoutFile.write("\n");
		else:
			eoutFile.write(fam);

	if etype != "mcl":
		eoutFile.write("\n##################################\n");
		if em == 1:
			eoutFile.write("# Rapid Contractions in ");
			eoutFile.write(specofint)
			eoutFile.write(" (");
			eoutFile.write(str(numc));
			eoutFile.write(")");
		elif em == 2:
			eoutFile.write("# ALL Contractions in ");
			eoutFile.write(specofint)
			eoutFile.write(" (");
			eoutFile.write(str(numc));
			eoutFile.write(")");
		eoutFile.write("\n##################################\n");

	if xtra != "":
		eoutFile.write(xtra);

	for fam in contracts:
		if etype == "mcl":
			fam = fam.replace("\n","").split("\t");
			for x in xrange(len(fam)):
				if x == 0:
					if fam[1] == '':
						famsize = "0";
					else:
						famsize = str(len(fam)-1);
					outline = fam[x] + " (" + famsize + " genes) -\n";
					eoutFile.write(outline);
				else:
					outline = fam[x] + "\n";
					eoutFile.write(outline);
			eoutFile.write("\n");
		else:
			eoutFile.write(fam);

	if etype == "report":
		if em == 1:
			print specofint + ": Rapidly changing families written in CAFE report format to " + Out;
		elif em == 2:
			print specofint + ": ALL changing families written in CAFE report format to " + Out;
	if etype == "mcl":
		if em == 1:
			print specofint + ": Rapidly changing families written in MCL dump format to " + Out;
		elif em == 2:
			print specofint + ": ALL changing families written in MCL dump format to " + Out;
	if etype == "cafe":
		if em == 1:
			print specofint + ": Rapidly changing families written in CAFE input format to " + Out;
		elif em == 2:
			print specofint + ": ALL changing families written in CAFE input format to " + Out;

############################################
#Main block
############################################

inFilename, soi, emode, outFilename, mclInFilename, mclOutFilename, cafeInFilename, cafeOutFilename, omode = optParse(0);
#Get the input parameters. soi is species of interest.

print "=======================================================================";
print "Analyzing file:\t\t\t\t" + inFilename;
if soi != "":
	print "Species of interest:\t\t\t" + soi;

	if omode == 0:
		print "-p 0\t\t\t\t\tWriting output for all species."
	else:
		print "-p 1\t\t\t\t\tWriting output for only species of interest.";

	outline = "Extract mode = " + str(emode) + ":\t\t\t";
	if emode == 1:
		outline = outline + "Extracting only families in which the species of interest is rapidly changing.";
	elif emode == 2:
		outline = outline + "Extracting all families in which the species of interest is changing.";
	print outline;

	if outFilename != None:
		print "Writing families in Report format to:\t" + outFilename;
	if mclInFilename != None and mclOutFilename != None:
		print "Using MCL input file:\t\t\t" + mclInFilename;
		print "Writing families in MCL format to:\t" + mclOutFilename;
	if cafeInFilename != None and cafeOutFilename != None:
		print "Using CAFE input file:\t\t\t" + cafeInFilename;
		print "Writing families in CAFE format to:\t" + cafeOutFilename;

	rapfilename = soi + "_rapids_ancestral.txt";
	print "Writing ancestral changes for\nspecies of interest to:\t\t\t" + rapfilename;

	if outFilename == None and mclInFilename == None and mclOutFilename == None and cafeInFilename == None and cafeOutFilename == None:
		print "------";
		print "Writing no output. Simply printing the analysis table for all species.";
else:
	print "------";
	print "Writing no output. Simply printing the analysis table for all species.";
print "=======================================================================";

############################################
#Parsing the CAFE report file
############################################

print "Reading and parsing CAFE report file...";

inFile = open(inFilename, "r");
inLines = inFile.readlines();
inFile.close();
#Reads the file to inLines

if inLines[2].find("Lambda tree:") != -1:
	treeline = inLines[0];
	nodeIDline = inLines[3];
	branchIDline = inLines[4];
	avgline = inLines[6];
	linestart = 10;
else:
	treeline = inLines[0];
	nodeIDline = inLines[2];
	branchIDline = inLines[3];
	avgline = inLines[5];
	linestart = 9;
#This checks if a Lambda tree line is present in the report file and adjusts how it gathers the tree,
#node IDs, and sister branch IDs appropriately.
#treeline: The original tree as input into CAFE.
#nodeIDline: The line from the report file which contains the node IDs used by CAFE.
#branchIDline: The line from the report file which pairs sister branches together.

avgline = avgline[avgline.index("\t")+1:];
avgline = avgline.replace("\n", "");
avgline = avgline.replace(")", "");
avgline = avgline.replace("(", "");
avgline = avgline.split("\t");
for b in xrange(len(avgline)):
	avgline[b] = avgline[b].split(",");
#Parsing of the average expansion line to be added to the results later.

branchIDline = branchIDline[branchIDline.index("=")+2:];
branchIDline = branchIDline[branchIDline.index(":")+2:];
branchIDline = branchIDline.replace("\n", "");
sisterPairs = branchIDline.split(" ");
sisterPairs = filter(None, sisterPairs);

sisterDict = {};

for b in xrange(len(sisterPairs)):
	siskey = sisterPairs[b][1:sisterPairs[b].index(",")];
	sisval = sisterPairs[b][sisterPairs[b].index(",")+1:sisterPairs[b].index(")")];
	pos = str(b+1) + ",1";
	sisterDict[siskey] = [];
	sisterDict[siskey].append(sisval);
	sisterDict[siskey].append(pos);

	siskey = sisterPairs[b][sisterPairs[b].index(",")+1:sisterPairs[b].index(")")];
	sisval = sisterPairs[b][1:sisterPairs[b].index(",")];
	pos = str(b+1) + ",2";
	sisterDict[siskey] = [];
	sisterDict[siskey].append(sisval);
	sisterDict[siskey].append(pos);
#The above block formats the branchIDline into sister pairs and stores the pairs in a dictionary of
#the key:value format node ID:[sister node ID, position in sister pair list].

treeline = treeline[5:]
treeline = treeline.replace("\n","");
nodeIDline = nodeIDline[15:];
nodeIDline = nodeIDline.replace("\n", "");
orignodeline = nodeIDline;
#Parsing of the tree and nodeIDlines. An original copy of the nodeIDline is kept for later reference.

root = nodeIDline[nodeIDline.rindex("<")+1:nodeIDline.rindex(">")];
#Identification of the root node label done by reverse searching the node ID tree (the root always appears
#at the end of the tree.

oldtree = treeline;
newtree = "";

for b in xrange(len(orignodeline)):

	if orignodeline[b] == ",":
		newtree = newtree + oldtree[oldtree.index(":"):oldtree.index(",")];
		oldtree = oldtree[oldtree.index(",")+1:];

	elif orignodeline[b] == ")":
		newtree = newtree + oldtree[oldtree.index(":"):oldtree.index(")")];
		oldtree = oldtree[oldtree.index(")")+1:];

	newtree = newtree + orignodeline[b];
newtree_copy = newtree;
#This block builds a newtree which incorporates the branch lengths from the original tree AND the node
#ID labels used by CAFE.

#print newtree;
#sys.exit();

nodeIDline = nodeIDline.replace(")", "&");
nodeIDline = nodeIDline.replace("(", "&");
nodeIDline = nodeIDline.replace(",", "&");
nodeIDline = nodeIDline.split("&");
nodeIDline = filter(None, nodeIDline);
#Additional parsing of the node ID line

tracker = {};
#The tracker dictionary is very important! It gathers all relevant information for each node in the following
#key:value format -- node ID:[branch length of node ID, ancestral node ID, branch length of ancestral node ID, position in sister pair list]
#If the node is a tip, the node ID is the species label.
#The position in the sister pair list is the same as those gathered above in the sisterPairs dictionary. It is a list of length 2 in the following
#format: [position of sister pair in list, position of node ID within sisterpair. ie -- if the sister pair list appears as: (0,1) (2,3), then the
#position in the sister pair list of node 2 is (2,1).

for each in nodeIDline:
	if len(each)-1 == each.index(">") - each.index("<"):
		key = each[each.index("<")+1:each.index(">")];
		if key != root:
			tracker[key] = [];

	else:
		key = each[:each.index("<")];
		tracker[key] = [];
#This loop builds the keys for the tracker based on all nodes in the node ID line.

rootnodes = [];
#rootbranches will store the branch lengths or the tip labels of the branches coming directly off of the root.

for key in tracker:
#This loop gathers the information in the tracker.

	if key.isdigit():
		curnode = key;
		curid = "<"+key+">";

	else:
		curnode = orignodeline[orignodeline.index(key)+len(key)+1:orignodeline.index(">", orignodeline.index(key))];
		curid = key;

	##Good for debugging
	##print "------";
	##print curid;
	##print "node", curnode;

	curbranch = getBranchLength(curid, newtree);
	##print "branch", curbranch;
	#Gets the branch length of the current node.

	tracker[key].append(curbranch);

	cursis = sisterDict[curnode][0];
	##print "sis", cursis;
	#Gets the sister node to find the ancestral node.

	if int(cursis) > int(curnode):
		anccheck = cursis;
	else:
		anccheck = curnode;

	anccheck = "<" + anccheck + ">"
	##print "check", anccheck;

	start = orignodeline.index(anccheck) + len(anccheck);

	curanc = orignodeline[orignodeline.index("<", start)+1:orignodeline.index(">", start)];
	ancid = "<" + curanc + ">";
	##print "anc", curanc;
	#Gets the ancestral node. The ancestral node is always direclty to the right of the larger node in each sister pair.

	tracker[key].append(curanc);

	if curanc == root:
		rootnodes.append(key);
		tracker[key].append("");
	else:
		ancbranch = getBranchLength(ancid, newtree);
		tracker[key].append(ancbranch);
	#If the ancestor is the root, there is no branch length associated with it, otherwise this gets the branch length
	#of the ancestral node. Also, if the ancestor is the root, this ID is added to the rootbranches list.

	tracker[key].append(sisterDict[curnode][1]);
	#Retrieves the position in the sister pair list of the current node.

##print tracker["gibbon"];
##sys.exit();

results = {};
#results is another important dictionary! It stores the counts of all family/gene gains and losses for each species. It
#has the following key:value format: species ID:[# fams expand, # genes expand, # fams equal, # fams contract, # genes contract, average expansion, # fams sigexpand, # fams sigcontract, total sig changes]

for tip in tracker:
	if tip.isdigit():
		continue;

	results[tip] = [0,0,0,0,0,0,0,0,0,0];
#Results are only tracked for the tips (species) of the tree. This initializes the dictionary accordingly.

for b in xrange(len(rootnodes)):
	if rootnodes[b].isdigit():
		rootnodes[b] = "<" + rootnodes[b] + ">";

#rootbranches are parsed here.

############################################
#Example tracker and results formats.
#tracker = {
#"Elephant":["101.7", "1", "", "1,1"],
#"Human":["94.2", "3", "7.5", "2,1"],
#"3":["7.5", "1", "", "1,2"],
#"Horse":["82.4", "5", "11.8", "3,1"],
#"5":["11.8", "3", "7.5", "2,2"],
#"Pig":["63.1", "7", "14.3", "4,1"],
#"7":["14.3", "9", "5", "5,1"],
#"Cow":["63.1", "7", "14.3", "4,2"],
#"9":["5", "5", "11.8", "3,2"],
#"Cat":["55.1", "11", "22.3", "6,1"],
#"11":["22.3", "9", "5", "5,2"],
#"Dog":["42.6", "13", "12.5", "7,1"],
#"13":["12.5", "11", "22.3", "6,2"],
#"Ferret":["38", "15", "4.6", "8,1"],
#"15":["4.6", "13", "12.5", "7,2"],
#"Panda":["38", "15", "4.6", "8,2"],
#}
#Node:[branch length,ancestral node, ancestral branch length, coordinates in cafe output]

#results = {
#"Elephant":[0,0,0,0,0,0,0,0],
#"Human":[0,0,0,0,0,0,0,0],
#"Horse":[0,0,0,0,0,0,0,0],
#"Pig":[0,0,0,0,0,0,0,0],
#"Cow":[0,0,0,0,0,0,0,0],
#"Cat":[0,0,0,0,0,0,0,0],
#"Dog":[0,0,0,0,0,0,0,0],
#"Ferret":[0,0,0,0,0,0,0,0],
#"Panda":[0,0,0,0,0,0,0,0]
#}
#[expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigequal,sigcontract,total sig changes]
############################################

############################################
#Thus ends the parsing of the input information. Now we can analyze the actual results...
############################################

if soi not in tracker:
	core.errorOut(7, "Species label for species of interest (-s) not found in input file. Re-check your labels!");
	optParse(1);

print "Analyzing CAFE report...";

numbars = 0;
donepercent = [];
#For the loading bar

famids = [];
soilosses = 0;
soilosses = [];
soiexpands = [];
soicontracts = [];

numrapids = 0;
soirapids = {};

#Initialization of species of interest counter variables.

for x in xrange(len(inLines)):

	numbars, donepercent = core.loadingBar(x, len(inLines), donepercent, numbars);
	#The loading bar.

	incount = 0;
	soicount = 1;
	outcount = 0;
	#Species of interest and its associated in and outgroups each have counters.

	if x > linestart:
	#Only get the lines which contain output for a family, therefore each line represents one family.
		line = inLines[x].split("\t");
		curfam = int(line[0]);
		curtree = line[1];
		curtree = curtree.replace("\n", "");
		newcurtree = ""
		famids.append(curfam);

		if float(line[2]) < 0.01:
			numrapids = numrapids + 1;

		for c in xrange(len(newtree)):
			newcurtree = newcurtree + newtree[c];
			if newtree[c] == ">":
				if c != len(newtree)-1:
					for b in xrange(len(curtree)):
						if curtree[b] == "_":
							curcount = curtree[b:curtree.index(":", b)];
							curtree = curtree[curtree.index(":", b)+1:];
							newcurtree = newcurtree + curcount;
							break;
		newcurtree = newcurtree + curtree[curtree.index("_"):];
		#This block builds a new tree out of the current output tree which incorporates both node IDs and
		#gene counts for each species.

		for each in tracker:
			if each.isdigit():
				curspec = "<" + each + ">";
			else:
				curspec = each;

			speccount = newcurtree[newcurtree.index("_", newcurtree.index(curspec))+1:newcurtree.index(":", newcurtree.index(curspec))];
			##print curspec;
#			if curspec == soi:
#			print speccount;
			#Getting the current species (or node) and its gene count for this family.


			curanc = "<" + tracker[each][1] + ">";

			if curspec in rootnodes:
				curanc = "root";

			if curanc != "root":
				anccount = newcurtree[newcurtree.index("_", newcurtree.index(curanc))+1:newcurtree.index(":", newcurtree.index(curanc))];
			else:
				anccount = newcurtree[newcurtree.rindex("_")+1:];
			#Then the ancestral gene count is retrieved to determine whether or not the current species has gained or lost genes in this family

			if not each.isdigit():
			#For the tips
				if anccount != "0" and speccount == "0":
					results[each][5] = results[each][5] + 1;

				if int(speccount) > int(anccount):
					curaction = "expand";
					results[each][0] = results[each][0] + 1;

					results[each][1] = results[each][1] + (int(speccount) - int(anccount));

					if curspec == soi and emode == 2:
						soiexpands.append(inLines[x]);

				elif int(speccount) == int(anccount):
					curaction = "equal";
					results[each][2] = results[each][2] + 1;

				elif int(speccount) < int(anccount):
					curaction = "contract";
					results[each][3] = results[each][3] + 1;

					results[each][4] = results[each][4] + (int(anccount) - int(speccount));

					if curspec == soi and emode == 2:
						soicontracts.append(inLines[x]);
				#Each species gene count is checked against its ancestral gene count to determine the change. These numbers are then
				#added to the results dictionary appropriately.

				if float(line[2]) < 0.01:
				#If it is determined that this family has changed significantly, then we check if this species has changed significantly.

					curpout = line[3];
					curpout = curpout.split(",(");

					for y in xrange(len(curpout)):
						curpout[y] = curpout[y].replace("(", "");
						curpout[y] = curpout[y].replace(")", "");
						curpout[y] = curpout[y].split(",");

					curpos = tracker[each][3].split(",");
					pos1 = int(curpos[0]);
					pos2 = int(curpos[1]);
					#In which case we use the position in the sister pair list to check the p-values output by CAFE for this family. The
					#position of the current species is retrieved from tracker and then the current p-value list is checked at that position.

					specpval = curpout[pos1-1][pos2-1];

					if float(specpval) < 0.01:
					#If that p-value is significant, this species has evolved rapidly in this family.
						if curspec == soi:
							if curfam not in soirapids:
								soirapids[curfam] = [];
							soirapids[curfam].append(anccount);
							soirapids[curfam].append(speccount);

						if curaction == "expand":
							results[each][7] = results[each][7] + 1;

							if curspec == soi and emode == 1:
								soiexpands.append(inLines[x]);

							if curspec == soi:
								soirapids[curfam].append("+");

						elif curaction == "contract":
							results[each][8] = results[each][8] + 1;

							if curspec == soi and emode == 1:
								soicontracts.append(inLines[x]);

							if curspec == soi:
								soirapids[curfam].append("-");
						#The significant change is categorized and the counts are added to results.


for key in results:

	curpos = tracker[key][3].split(",");
	pos1 = int(curpos[0])-1;
	pos2 = int(curpos[1])-1;
	results[key][6] = avgline[pos1][pos2]

	results[key][9] = results[key][8] + results[key][7];
#Getting the average expansion for each species and the total number of significant changes.

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print "=======================================================================";
#Loading bar

############################################
#Output
############################################

print "\tExpansions\tGenes Gained\tEqual\tContractions\tGenes Lost\tFamilies Lost\tAverage Expansion\tSig Expansions\tSig Contractions\tTotal Sig Changes";
for species in results:
	if omode == 1 and species != soi:
		continue;
	outline = species + "\t";
	for col in results[species]:
		outline = outline + str(col) + "\t";
	print outline;
#This block simply prints the information stored in results to the screen.

if soi != 0:
#Output for the species of interest, if defined.
	print "------";

	rFile = open(rapfilename, "w");
	outline = "CAFE family ID\tType of change\tAncestral gene count\t" + soi + " gene count\tDifference\n";
	rFile.write(outline);

	for each in famids:
		if each in soirapids:
			curdiff = soirapids[each][2] + str(abs(int(soirapids[each][0]) - int(soirapids[each][1])));
			outline = str(each) + "\t" + soirapids[each][2] + "\t" + soirapids[each][0] + "\t" + soirapids[each][1] + "\t" + curdiff + "\n";
			rFile.write(outline);
	rFile.close();

	if emode != 0:
	#If an extract mode is specified. (Remember, -e 1 extracts only rapidly changing families for the species of interest, -e 2 extracts
	#all changing families from soi.
		if emode == 1:
			numsoiexp = results[soi][6];
			numsoicon = results[soi][7];
		elif emode == 2:
			numsoiexp = results[soi][0];
			numsoicon = results[soi][3];

		if outFilename != None:

			extra = "";
			for line in inLines:
				if line[:7] == "Average":
					break;
				extra = extra + line;

			e_type = "report";
			writeRapids(inFilename, outFilename, soiexpands, soicontracts, soi, inFilename, emode, numsoiexp, numsoicon, e_type, extra);
		#This block calls the write function if the families are extracted in CAFE report format (if -o is defined).

		if mclInFilename != None:

			mclexpands = [];
			mclcontracts = [];
			extra = ""

			mclFile = open(mclInFilename, "r");
			mclLines = mclFile.readlines();
			mclFile.close();

			for fam in soiexpands:

				fam = fam.split("\t");
				fam_index = int(fam[0])-1;

				efam = mclLines[fam_index];

				if omode == 1:
					newefam = [];
					efam = efam.split("\t");
					for each in efam:
						if each.find(soi) != -1:
							newefam.append(each);
					efam = "\t".join(newefam) + "\n";
					if efam.find("\n") == -1:
						efam = efam + "\n";

				outline = str(fam_index+1) + "\t" + efam;
				mclexpands.append(outline);

			for fam in soicontracts:

				fam = fam.split("\t");
				fam_index = int(fam[0])-1;

				mclfam = mclLines[fam_index];

				if omode == 1:
					newmclfam = [];
					mclfam = mclfam.split("\t");
					for each in mclfam:
						if each.find(soi) != -1:
							newmclfam.append(each);
					mclfam = "\t".join(newmclfam);
					if mclfam.find("\n") == -1:
						mclfam = mclfam + "\n";

				outline = str(fam_index+1) + "\t" + mclfam;
				mclcontracts.append(outline);
			del mclLines;

			e_type = "mcl";
			writeRapids(mclInFilename, mclOutFilename, mclexpands, mclcontracts, soi, inFilename, emode, numsoiexp, numsoicon, e_type, extra);
		#This block extracts the families from an MCL dump file and passes them to the write function (if -m and -u are defined).

		if cafeInFilename != None:
			rotator = 0;
			cinexpands = [];
			cincontracts = [];
			extra = "";

			cinwork = "Extracting families from CAFE input file...  ";
			sys.stderr.write(cinwork);

			cinFile = open(cafeInFilename, "r");
			cinLines = cinFile.readlines();
			cinFile.close();
			extra = cinLines[0];

			for fam in soiexpands:

				fam = fam.split("\t");
				fam_index = fam[0];

				for b in xrange(len(cinLines)):
					rotator = core.loadingRotator(b, rotator, 1000);
					if b == 0:
						continue;

					cinfam = cinLines[b].split("\t");
					if cinfam[1] == fam_index:
						cinfam = "\t".join(cinfam);
						cinexpands.append(cinfam);
						break;

			for fam in soicontracts:

				fam = fam.split("\t");
				fam_index = fam[0];

				for b in xrange(len(cinLines)):
					rotator = core.loadingRotator(b, rotator, 1000);
					if b == 0:
						continue;

					cinfam = cinLines[b].split("\t");
					if cinfam[1] == fam_index:
						cinfam = "\t".join(cinfam);
						cincontracts.append(cinfam);
						break;
			del cinLines;

			sys.stderr.write('\b' * (len(cinwork) + 1));

			e_type = "cafe";
			writeRapids(cafeInFilename, cafeOutFilename, cinexpands, cincontracts, soi, inFilename, emode, numsoiexp, numsoicon, e_type, extra);
		#This block extracts the families from the CAFE input file and passes them to the write function (if -c and -t are defined).
print "The total number of rapidly evolving families:\t", numrapids;
print "=======================================================================";
