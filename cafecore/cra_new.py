#!/usr/bin/python
#############################################################################
# CRA, or CAFE Report Analysis, takes as input a report file from a CAFE run
# and displays the results in a more interpretable manner. A table of counts
# of gene family data is printed out for each node and an output file is specified
# that lists rapidly evolving families across the tree and for each node.
#
#Usage: python cra_new.py -i [a CAFE report file (.cafe)] -o [output file name]
#
#Depenencies: core, treeparse
#
#Gregg Thomas, Spring 2016
#############################################################################

import sys, os, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../corelib/"))
import core, treeparse

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="report_input_file", help="A CAFE report file (.cafe).");
	parser.add_argument("-o", dest="output_file", help="An output file to write a list of all rapidly evolving families for each node.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.report_input_file == None or args.output_file == None:
			core.errorOut(1, "A CAFE report file must be specified with -i");
			optParse(1);

		return args.report_input_file, args.output_file;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

#######################

def formatLineParse(line):
# This function handles CAFE's weird node pair format for the p-values and node ids.

	if "=" in line:
		line = line.split("=")[1];
	if ":" in line:
		line = line.split(": ")[1].strip();
	line = line.replace("(", "").replace(")", "");
	line = line.split(" ");
	line = [f.split(",") for f in line];

	return line;

#######################

def nodeRelabel(treedict):
# Family trees are read with gene counts on the tip labels. This function removes them.

	tmp = {};
	#print treedict;

	for oldkey in treedict:
		if treedict[oldkey][3] == 'tip':
			newkey = oldkey[:oldkey.index("_")];
			tmp[newkey] = treedict[oldkey];
		else:
			tmp[oldkey] = treedict[oldkey];

	return tmp;

#######################

def nodeMap(cafetd, mytd):
# CAFE has pre-determined mappings in the tree. When I read the tree with my own script the mappings
# are different. This function creates a map from my node ids to CAFE's node ids.
# The dictionary nodemap has the following {key:value} format: {my node id:CAFE's node id}

	nodemap = {};
	# The map dictionary.

	##############
	# for node in cafetd:
	# 	if cafetd[node][3] == 'tip':
	# 		spec = node[:node.index("<")];
	# 		cafeid = node[node.index("<")+1:node.index(">")];
	# 		nodemap[cafeid] = spec;

	# while len(nodemap) != len(cafetd):
	# 	for node in cafetd:
	# 		if cafetd[node][3] == 'root':
	# 			continue;

	# 		orignode = node;
	# 		node = node[node.index("<")+1:node.index(">")];

	# 		if node in nodemap:
	# 			if cafetd[orignode][3] == 'tip':
	# 				curanc = cafetd[orignode][1];
	# 				mapanc = mytd[nodemap[node]][1];
	# 			else:
	# 				curanc = cafetd[orignode][1];
	# 				mapanc = mytd["<" + nodemap[node] + ">"][1];

	# 			nodemap[curanc.replace("<","").replace(">","")] = mapanc.replace("<","").replace(">","");
	##############
	# The above formats nodemap with the reverse {key:value} format: {CAFE's node id:my node id}

	for node in cafetd:
		if cafetd[node][3] == 'tip':
			spec = node[:node.index("<")];
			cafeid = node[node.index("<"):];
			nodemap[spec] = node;
	# First map the tips by their unique species labels.

	while len(nodemap) != len(mytd):
		for node in mytd:
			if mytd[node][3] == 'root':
				continue;

			if node in nodemap:
				curanc = mytd[node][1];
				mapanc = cafetd[nodemap[node]][1];

				nodemap[curanc] = mapanc;
	# Then do a post-order traversal and map the current node's ancestor to it's map's ancestor.

	return nodemap;


############################################
#Main block
############################################

infilename, outfilename = optParse(0);
#Get the input parameters.

print "=======================================================================";
print "\t\tCAFE Report File Analysis"
print "\t\t" + core.getDateTime();
print "---------";
print "Parsing format information...\n";

infile = open(infilename, "r");
inlines = infile.readlines();
infile.close();
# Reads the input report file.

if inlines[2].find("Lambda tree:") != -1:
	treeline = inlines[3];
	formatline = inlines[4];
	avgline = inlines[6];
	linestart = 11;
else:
	treeline = inlines[2]
	formatline = inlines[3];
	avgline = inlines[5];
	linestart = 10;
# If CAFE was run with a lambda tree structure as input, the report file places that on the third
# line. This shifts all the other relevant lines down by 1. This if/else accounts for that.

labeled_tree = treeline[treeline.index(":")+1:].strip();
tinfo, newtree = treeparse.treeParseNew(labeled_tree,2);
# This reads the CAFE tree with its node labels.

formatline = formatLineParse(formatline);
# formatline is CAFE's line with its paired node format with node ids. The formatLineParse function
# reads that format and returns it as a list of lists.

avgline = avgline.split(":\t")[1].strip().replace("\t", " ");
avgline = formatLineParse(avgline);
# The line of average expansions for each node, in the paired node format. Again passed to formatLineParse
# to make it interpretable.

print "---------";
print "Initializing output structures...\n";

outfile = open(outfilename, "w");
outfile.write("");
# Initialize the output file.

rapids = {"total" : []};
results = {};
for node in tinfo:
	if tinfo[node][3] == 'root':
		continue;
	rapids[node] = [[],[]];
	results[node] = [0,0,0,0,0,0,0,0,0,0];
	# [expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]
# rapids and results are the two main dictionaries to store CAFE's results.
# rapids {key:value} format: {node:list of two lists containing family ids for rapid expansions and rapid contractions, respectively}
# results {key:value} format: {node:[expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]}

for j in range(len(formatline)):
	for k in range(len(formatline[j])):
		n = "<" + formatline[j][k] + ">";
		for r in results:
			if n in r:
				results[r][6] = avgline[j][k];
# Setting average expansion in results as read from avgline

print "---------";
print "Counting changes per branch...\n";

numbars = 0;
donepercent = [];
i = 0;
acount = 0;

for inline in inlines:
# Each line of the report file is read.
	numbars, donepercent = core.loadingBar(i, len(inlines), donepercent, numbars);
	i = i + 1;

	if i <= linestart:
		continue;
	# If the line is a CAFE info line, skip it.

	inline = inline.strip().split("\t");
	famid = inline[0];
	famtree = inline[1];
	nodeformat = inline[3].replace("),(", ") (");
	# Parsing the information for the current family.

	tlinfo, newfamtree = treeparse.treeParseNew(famtree, 1);
	# Reading the tree and adding my node labels.

	for tlnode in tlinfo:
		if tlinfo[tlnode][3] == 'root':
			tlinfo[tlnode].append(famtree[famtree.rfind("_")+1:]);
		elif tlinfo[tlnode][3] == 'tip':
			tlinfo[tlnode].append(tlnode[tlnode.index("_")+1:]);
		else:
			tlinfo[tlnode][4] = tlinfo[tlnode][4][1:];
	# Gene counts for each node are read as support values for internal nodes, but must
	# have the underscore removed. Tip and root node counts are added here as well.

	tlinfo = nodeRelabel(tlinfo);
	# Removes the gene counts from the tip node labels.

	if i == (linestart + 1):
		maps = nodeMap(tinfo, tlinfo);
	# If this is the first family, we need to build our maps from my node ids to CAFE's.

	for tlnode in tlinfo:
	# For each node in the current gene family tree, we make our counts.

		if tlinfo[tlnode][3] == 'root':
			continue;
		# No count is made at the root of the tree.

		curanc = tlinfo[tlnode][1];
		curmap = maps[tlnode];
		# Get the ancestor and the map of the current node.

		curcount = int(tlinfo[tlnode][4]);
		anccount = int(tlinfo[curanc][4]);
		# Get the gene counts of the current node and the ancestor.

		diff = curcount - anccount;
		# Calculate the difference in gene count.

		typeflag = 0;
		# typeflag tells whether an expansion or contraction has occurred.

		if curcount > anccount:
			typeflag = 1;
			results[curmap][0] += 1;
			results[curmap][1] += diff;
		# If the difference in gene count between the current node and the ancestor is positive, an
		# expansion has occurred. This makes the appropriate counts.

		elif curcount < anccount:
			typeflag = 2
			results[curmap][3] += 1;
			results[curmap][4] += abs(diff);

			if curcount == 0 and anccount != 0:
				results[curmap][5] += 1;
		# If the difference in gene count between the current node and the ancestor is negative, a
		# contraction has occurred. This makes the appropriate counts. It also checks for family losses
		# along that branch by seeing if the current node has transitioned to a count of 0.

		elif curcount == anccount:
			results[curmap][2] += 1;
		# Otherwise, the counts at the current node and the ancestor are the same and no change has occurred.

		if float(inline[2]) < 0.01:
		# If the family p-value is below a threshold, the family is rapidly evolving.

			if famid not in rapids['total']:
				rapids['total'].append(famid);
				nodes = formatLineParse(nodeformat);
			# Add the family id to the 'total' key of rapids. This also parses the nodeformat line which is
			# in the paired CAFE node format.
			
			pairnodeid = curmap[curmap.index("<")+1:curmap.index(">")];
			# Since the paired format does not include the brackets that the other node labels do, I have to
			# remove them to check against that format.

			for j in range(len(nodes)):
				for k in range(len(nodes[j])):
					if formatline[j][k] == pairnodeid and float(nodes[j][k]) < 0.01:
						if typeflag == 1:
							results[curmap][7] += 1;
							rapids[curmap][0].append(famid);
						elif typeflag == 2:
							results[curmap][8] += 1;
							rapids[curmap][1].append(famid);
						results[curmap][9] += 1;

			# Runs through the paired format as a list of lists. If the p-value of that node is less than a threshold
			# that branch is rapidly evolving. Based on typeflag, the appropriate counts are made. The family id is
			# also added to the current node in rapids.

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print "=======================================================================";

outfile.write("# The labeled CAFE tree:\t" + labeled_tree + "\n");

outline = "Overall rapids:\t"
for f in rapids['total']:
	outline = outline + f + ",";
outline = outline[:-1] + "\n";
outfile.write(outline);

for spec in rapids:
	if spec == 'total':
		continue;

	for f in range(len(rapids[spec])):
		if f == 0:
			outline = spec + " rapid expansions:\t";
		elif f == 1:
			outline = spec + " rapid contractions:\t";

		for rapid_f in rapids[spec][f]:
			outline = outline + rapid_f + ",";
		outline = outline[:-1] + "\n";
		outfile.write(outline);

print "RESULTS TABLE -- tab delimted for easy copy/pasting into your favorite spreadsheet program"
print "\tExpansions\tGenes Gained\tEqual\tContractions\tGenes Lost\tFamilies Lost\tAverage Expansion\tSig Expansions\tSig Contractions\tTotal Sig Changes";
for species in results:
	outline = species + "\t";
	for col in results[species]:
		outline = outline + str(col) + "\t";
	print outline;

print
print "Number of rapidly evolving families in tree:\t", len(rapids['total']);
print "CAFE labeled tree:\t" + labeled_tree;
# This block simply prints the information stored in results to the screen.


