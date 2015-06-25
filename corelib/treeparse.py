#!/usr/bin/python
#############################################################################
#The treeParse function takes as input a phylogenetic tree with branch lengths
#and returns a dictionary with usable info about the tree in the following
#format:
#
#node:[branch length, ancestral node, ancestral branch length, sister node, sister branch length, descendent 1, descendent 1 branch length, descendent 2, descendent 2 branch length, node type]
#
#############################################################################

import sys
'''
import argparse

#############################################################################

def IO():

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_data", help="Either a tree or a file containing a tree.");
	parser.add_argument("-t", dest="input_type", help="Either 't' for tree or 'f' for file containing a tree.");
	parser.add_argument("-o", dest="output_file", help="Output file name.");

	args = parser.parse_args();

	if args.input_data == None:
		parser.print_help();
		sys.exit();

	if args.input_type not in ['t','f']:
		print " -------------------------------------------------------------------------------------------";
		print "|**Error 1: -t must take values of either 't' or 'f' depending on input type (Tree or File) |";
		print " -------------------------------------------------------------------------------------------";
		parser.print_help();
		sys.exit();

	return args.input_data, args.input_type, args.output_file;
'''
#############################################################################
def getBranchLength(bltree, spec_label):

	#print spec_label;
	d = 0;
	startind = 0;

	while d < (len(bltree)-1):
		if bltree[d] == ":":
			current_node = bltree[max(bltree.rfind("(",startind,d),bltree.rfind(")",startind,d),bltree.rfind(",",startind,d))+1:d];
			#print current_node;
			if current_node == spec_label:

				opind = bltree.find("(",d);
				cpind = bltree.find(")",d);
				coind = bltree.find(",",d);

				indcheck = [opind,cpind,coind];
				#print indcheck;

				for a in xrange(len(indcheck)):
					if indcheck[a] == -1:
						indcheck[a] = 10000;

				#print indcheck

				curbranch = bltree[d+1:min(indcheck)];
				return curbranch;
		d = d + 1;
	startind = d;

#############################################################################
def treeParse(tree):

	tree = tree.replace("\n","");

	if tree[len(tree)-1:] != ";":
		tree = tree + ";";

	new_tree = "";
	z = 0;
	numnodes = 1;

	while z < (len(tree)-1):
		if tree[z] == ":" and tree[z-1] == ")":
			new_tree = new_tree + "<" + str(numnodes) + ">";
			numnodes = numnodes + 1;
		new_tree = new_tree + tree[z];
		z = z + 1;
	rootnode = "<" + str(numnodes) + ">"
	new_tree = new_tree + rootnode;


	#print new_tree;
	#print "-----------------------------------";

	ancs = {};
	nofo = {};

	z = 0;
	startind = 0;
	while z < (len(new_tree)-1):

		if new_tree[z] == ":":
			curnode = new_tree[max(new_tree.rfind("(",startind,z),new_tree.rfind(")",startind,z),new_tree.rfind(",",startind,z))+1:z];
			numcpneeded = 1
			numcp = 0;
			nofo[curnode] = [];

			a = z;

			while a < (len(new_tree)-1):
				if new_tree[a] == "(":
					numcpneeded = numcpneeded + 1;
				if new_tree[a] == ")" and numcpneeded != numcp:
					numcp = numcp + 1;
				if new_tree[a] == ")" and numcpneeded == numcp:
					if a == (len(new_tree)-4):
						curanc = new_tree[a+1:];
					elif new_tree[a+1:].find(":") == -1:
						curanc = new_tree[len(new_tree)-4:];
					else:
						curanc = new_tree[a+1:new_tree.index(":", a)];
					a = 10000;

					ancs[curnode] = curanc;
				a = a + 1;

			startind = z;
		z = z + 1;

	#for key in ancs:
	#	print key + ":", ancs[key]
	#print "---------";

	nofo[rootnode] = [];

	for node in nofo:
		cur_bl = getBranchLength(new_tree,node);
		nofo[node].append(cur_bl);

		if node != rootnode:
			cur_anc = ancs[node];
			nofo[node].append(cur_anc);
			cur_anc_bl = getBranchLength(new_tree,cur_anc);
			nofo[node].append(cur_anc_bl);
			for each in ancs:
				if each != node and ancs[each] == cur_anc:
					cur_sis = each;
					nofo[node].append(cur_sis);
					cur_sis_bl = getBranchLength(new_tree,cur_sis);
					nofo[node].append(cur_sis_bl);
		else:
			nofo[node].append("");
			nofo[node].append("");
			nofo[node].append("");
			nofo[node].append("");

		tipflag = 1;

		for each in ancs:
			if ancs[each] == node:
				tipflag = 0;
				cur_desc = each;
				nofo[node].append(cur_desc);
				cur_desc_bl = getBranchLength(new_tree,cur_desc);
				nofo[node].append(cur_desc_bl);

		if tipflag == 1:
			nofo[node].append("");
			nofo[node].append("");
			nofo[node].append("");
			nofo[node].append("");

		if nofo[node][8] == "":
			nofo[node].append("tip");
		elif nofo[node][0] == None:
			nofo[node].append("root");
		else:
			nofo[node].append("internal");

	return nofo, new_tree;

#############################################################################
'''
#indataname, intype, outfilename = IO();

if intype == 't':
	intree = indataname;
elif intype == 'f':
	inFile = open(indataname, "r");
	intree = inFile.read();
	inFile.close();

node_info, ntree = treeParse(intree);

print "Node\tBranch length\tAncestral node\tAncestral branch length\tSister node\tSister branch length\tDescendent 1\t Descendent 1 branch length\tDescendent 2\tDescendent 2 branch length\tNode type";

for key in node_info:
	outline = key + "\t";
	for x in xrange(len(node_info[key])):
		if x != len(node_info[key])-1:
			outline = outline + str(node_info[key][x]) + "\t";
		else:
			outline = outline + str(node_info[key][x]);
	print outline;

if outfilename != None:
	outFile = open(outfilename, "w");
	outline = "Node\tBranch length\tAncestral node\tAncestral branch length\tSister node\tSister branch length\tDescendent 1\t Descendent 1 branch length\tDescendent 2\tDescendent 2 branch length\tNode type\n";
	outFile.write(outline);

	for key in node_info:
		outline = key + "\t";
		for x in xrange(len(node_info[key])):
			if x != len(node_info[key])-1:
				outline = outline + str(node_info[key][x]) + "\t";
			else:
				outline = outline + str(node_info[key][x]) + "\n";
		outFile.write(outline);
	outFile.close();

#infilename = "tree_parse_in.txt";

#inFile = open(infilename, "r");
#intree = inFile.read();
#inFile.close();

#node_info = treeParse(intree);

#outFile = open("testout.txt","w");
#outline = "Node\tBranch length\tAncestral node\tAncestral branch length\tSister node\tSister branch length\tDescendent 1\t Descendent 1 branch length\tDescendent 2\tDescendent 2 branch length\n";
#outFile.write(outline);

#for key in nofo:
#	outline = key + "\t";
#	for each in nofo[key]:
#		outline = outline + str(each) + "\t";
#	outFile.write(outline);
#	outFile.write("\n");
#outFile.close();
'''


