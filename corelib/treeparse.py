#!/usr/bin/python
#############################################################################
#Functions for retrieving information from newick formatted, rooted phylogenetic trees.
#Gregg Thomas
#Spring 2013-present
#############################################################################

import sys

#############################################################################
def getBranchLength(bltree, spec_label):
#Returns the branch length of a species given a newick formatted tree. Used by treeParse.

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

def comAnc(spec_list, treedict):
#Given a list of species within the tree and the dictionary returned by treeParse using that tree,
#this function checks whether those species are monophyletic (ie they all share a common ancestor).

	cur_list = [];
	for b in spec_list:
		if b in treedict:
			cur_list.append(b);

#	print treedict;
#	print spec_list;
#	print cur_list;

	ancdict = {};
	for b in cur_list:
		ancdict[b] = treedict[b][1];

	new_list = [];
	for b in ancdict:
		if ancdict.values().count(ancdict[b]) > 1 and ancdict[b] not in new_list:
			new_list.append(ancdict[b]);
		elif treedict[b][1] not in new_list:
			new_list.append(b);

	if not all(n in cur_list for n in new_list):
		flag = comAnc(new_list, treedict);
	elif len(new_list) > 1:
		#print "not monophyletic";
		flag = 0;
	else:
		#print "monophyletic";
		flag = 1;
		
	return flag;

#############################################################################

def treeParse(tree):
#The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a 
#dictionary with usable info about the tree in the following format:
#
#node:[branch length, ancestral node, ancestral branch length, sister node, sister branch length, descendent 1, descendent 1 branch length, descendent 2, descendent 2 branch length, node type]

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



