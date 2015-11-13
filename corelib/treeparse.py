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
			if current_node == spec_label:

				opind = bltree.find("(",d);
				cpind = bltree.find(")",d);
				coind = bltree.find(",",d);

				indcheck = [opind,cpind,coind];

				for a in xrange(len(indcheck)):
					if indcheck[a] == -1:
						indcheck[a] = 10000;

				curbranch = bltree[d+1:min(indcheck)];
				return curbranch;
		d = d + 1;
	startind = d;

#############################################################################

def getDesc(d_spec, d_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds the direct descendants of the species.
	d_list = [];
	for node in d_treedict:
		if d_treedict[node][1] == d_spec:
			d_list.append(node);

	if d_list == []:
		return [d_spec];
	else:
		return d_list;

#############################################################################

def getClade(c_spec, c_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds all tip labels that are descendants of the current node.
# This is done by getting the direct descendants of the current node with getDesc and then
# recursively calling itself on those descendants.
	clade = [];
	c_desc = getDesc(c_spec, c_treedict);
	for d in c_desc:
		if c_treedict[d][3] != 'tip':
			clade.append(getClade(d, c_treedict));
		else:
			clade.append(d);

	r_clade = [];
	for c in clade:
		if type(c) == list:
			for cc in c:
				r_clade.append(cc);
		else:
			r_clade.append(c);

	return r_clade;

#############################################################################

def LCA(spec_list, treedict):
	#print treedict;
	#print spec_list;
	#if spec_list.count(spec_list[0]) == len(spec_list):
	#	return treedict[spec_list[0]][1], 1;

	ancs = {};
	for spec in spec_list:
		ancs[spec] = [spec];

	for spec in spec_list:
		if treedict[spec][3] == 'root':
			continue;
		curanc = treedict[spec][1];
		ancs[spec].append(curanc);
		while treedict[curanc][3] != 'root':
			curanc = treedict[curanc][1];
			ancs[spec].append(curanc);
	#print ancs;

	intersect_anc = set.intersection(*map(set, ancs.values()))
	lcp = [t for t in ancs.values()[0] if t in intersect_anc]

	#lcp = sorted(set.intersection(*map(set, ancs.values())), key=lambda x: ancs.values()[0].index(x))
	monophyletic = 0;
	if set(getClade(lcp[0],treedict)) == set(spec_list):
		monophyletic = 1;

	#print ancs;
	#print ancs;
	#print getClade(lcp[0],treedict);
	#print monophyletic;
	#sys.exit();
	return lcp[0], monophyletic;

#############################################################################

def comAnc(spec_list, treedict):
#Given a list of species within the tree and the dictionary returned by treeParse using that tree,
#this function checks whether those species are monophyletic (ie they all share a common ancestor).

	#print spec_list;
	if len(spec_list) > 1:
		if spec_list.count(spec_list[0]) == len(spec_list):
			if treedict[spec_list[0]][3] == 'root':
				return 1, spec_list[0];
			else:
				return 1, treedict[spec_list[0]][1];
		if treedict[spec_list[0]][1] == spec_list[1]:
			return 1, spec_list[1];
		if treedict[spec_list[1]][1] == spec_list[0]:
			return 1, spec_list[0];

	cur_list = [];
	for b in spec_list:
		if b in treedict:
			cur_list.append(b);

	ancdict = {};
	for b in cur_list:
		ancdict[b] = treedict[b][1];

	new_list = [];
	for b in ancdict:
		if ancdict.values().count(ancdict[b]) > 1 and ancdict[b] not in new_list:
			new_list.append(ancdict[b]);
		elif treedict[b][1] not in new_list:
			new_list.append(b);

	#print new_list;

	if not all(n in cur_list for n in new_list):
		flag, com_anc = comAnc(new_list, treedict);
	elif len(new_list) > 1:
		#print "not monophyletic";
		flag = 0;
		com_anc = "";
	else:
		#print "monophyletic";
		flag = 1;
		com_anc = new_list[0];

	return flag, com_anc;

#############################################################################

def nodeDepth(n_spec, n_treedict):
	if n_treedict[n_spec][3] == 'root':
		return [];

	ancs = []
	curanc = n_treedict[n_spec][1];
	ancs.append(curanc);
	while n_treedict[curanc][3] != 'root':
		curanc = n_treedict[curanc][1];
		ancs.append(curanc);
	return ancs;

#############################################################################

def specRelabel(s, t_d):
#Relabels species to match the labels in the tree. (i5k)
	for node in t_d:
		if s in node:
			s = node;
	return s;

#############################################################################

def numInternal(treedict):
	num_nodes = 0;
	for node in treedict:
		if treedict[node][3] != 'tip':
			num_nodes = num_nodes + 1;
	return num_nodes;

#############################################################################

def treeParseNew(tree, tree_type):
#The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
#dictionary with usable info about the tree in the following format:
#New (current) format:
#node:[branch length, ancestral node, ancestral branch length, node type]
#
#Old format
#node:[branch length, ancestral node, ancestral branch length, sister node, sister branch length, descendent 1, descendent 1 branch length, descendent 2, descendent 2 branch length, node type]
#
#
#Tree type 1: tree has branch lengths.
#Tree type 2: tree is just topology.

	tree = tree.replace("\n","");
	if tree[len(tree)-1:] != ";":
		tree = tree + ";";
	##Some string handling

	new_tree = "";
	z = 0;
	numnodes = 1;

	while z < (len(tree)-1):
		if tree_type == 1:
			if tree[z] == ":" and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		if tree_type == 2:
			if (tree[z] == "," or tree[z] == ")") and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		new_tree = new_tree + tree[z];
		z = z + 1;
	if new_tree[-1] == ")":
		rootnode = "<" + str(numnodes) + ">"
		new_tree = new_tree + rootnode;
	else:
		rootnode = new_tree[new_tree.rfind(")")+1:];

	##This first block labels all internal nodes with the format <#>

	#print new_tree;
	#print "-----------------------------------";

	ancs = {};
	nofo = {};

	z = 0;
	startind = 0;
	while z < (len(new_tree)-1):
	##Here, the ancestral nodes of each node are found

		if tree_type == 1:
		##The major difference between trees with branch lengths (type 1) and without (type 2) is seen here. Finding the ancestral nodes requires
		##almost a completely different set of logic statements.
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
						a = 100000000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		if tree_type == 2:
			if new_tree[z] == "," or new_tree[z] == ")":
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
						else:
							mindex = 999999999;
							for c in ["(",")",","]:
								cind = new_tree.find(c,a+1);
								if cind < mindex and cind != -1:
									mindex = cind;
									minchar = c;
							curanc = new_tree[a+1:mindex];
						a = 10000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		z = z + 1;
	##End ancestral node block
#	print curanc;
	for key in ancs:
		print key + ":", ancs[key]
	print "---------";

	##The next block gets all the other info for each node: sister and decendent nodes and branch lengths (if type 1)
	##and node type (tip, internal, root). This is easy now that the ancestral nodes are stored.
	nofo[rootnode] = [];
	for node in nofo:
		if tree_type == 1:
			cur_bl = getBranchLength(new_tree,node);
		elif tree_type == 2:
			if node == rootnode:
				cur_bl = None;
			else:
				cur_bl = "NA";
		nofo[node].append(cur_bl);

		if node != rootnode:
			cur_anc = ancs[node];
			nofo[node].append(cur_anc);
			if tree_type == 1:
				cur_anc_bl = getBranchLength(new_tree,cur_anc);
			elif tree_type == 2:
				cur_anc_bl = "NA";
			nofo[node].append(cur_anc_bl);
		else:
			j = 0;
			while j < 2:
				nofo[node].append("");
				j = j + 1;

	for node in nofo:
		if node == rootnode:
			nofo[node].append("root");
		elif getDesc(node,nofo) == [node]:
			nofo[node].append("tip");
		else:
			nofo[node].append("internal");

	##End info retrieval block.

#	for key in nofo:
#		print key + ":" + str(nofo[key]);


	return nofo, new_tree;

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
# THIS IS THE OLD TREEPARSE FUNCTION that I am keeping around until I can make the switch to the
# new one in all my scripts.
def treeParse(tree, tree_type):
#The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
#dictionary with usable info about the tree in the following format:
#
#node:[branch length, ancestral node, ancestral branch length, sister node, sister branch length, descendent 1, descendent 1 branch length, descendent 2, descendent 2 branch length, node type]
#
#New proposed format:
#node:[branch length, ancestral node, ancestral branch length, node type]
#Tree type 1: tree has branch lengths.
#Tree type 2: tree is just topology.

	tree = tree.replace("\n","");
	if tree[len(tree)-1:] != ";":
		tree = tree + ";";
	##Some string handling

	new_tree = "";
	z = 0;
	numnodes = 1;

	while z < (len(tree)-1):
		if tree_type == 1:
			if tree[z] == ":" and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		if tree_type == 2:
			if (tree[z] == "," or tree[z] == ")") and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		new_tree = new_tree + tree[z];
		z = z + 1;
	if new_tree[-1] == ")":
		rootnode = "<" + str(numnodes) + ">"
		new_tree = new_tree + rootnode;
	else:
		rootnode = new_tree[new_tree.rfind(")")+1:];

	##This first block labels all internal nodes with the format <#>

#	print new_tree;
#	print "-----------------------------------";

	ancs = {};
	nofo = {};

	z = 0;
	startind = 0;
	while z < (len(new_tree)-1):
	##Here, the ancestral nodes of each node are found
		if tree_type == 1:
		##The major difference between trees with branch lengths (type 1) and without (type 2) is seen here. Finding the ancestral nodes requires
		##almost a completely different set of logic statements.
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

		if tree_type == 2:
			if new_tree[z] == "," or new_tree[z] == ")":
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
						else:
							mindex = 999999999;
							for c in ["(",")",","]:
								cind = new_tree.find(c,a+1);
								if cind < mindex and cind != -1:
									mindex = cind;
									minchar = c;
							curanc = new_tree[a+1:mindex];
						a = 10000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		z = z + 1;
	##End ancestral node block
#	print curanc;
#	for key in ancs:
#		print key + ":", ancs[key]
#	print "---------";

	##The next block gets all the other info for each node: sister and decendent nodes and branch lengths (if type 1)
	##and node type (tip, internal, root). This is easy now that the ancestral nodes are stored.
	nofo[rootnode] = [];

	for node in nofo:
		if tree_type == 1:
			cur_bl = getBranchLength(new_tree,node);
		elif tree_type == 2:
			if node == rootnode:
				cur_bl = None;
			else:
				cur_bl = "NA";
		nofo[node].append(cur_bl);

		if node != rootnode:
			cur_anc = ancs[node];
			nofo[node].append(cur_anc);
			if tree_type == 1:
				cur_anc_bl = getBranchLength(new_tree,cur_anc);
			elif tree_type == 2:
				cur_anc_bl = "NA";
			nofo[node].append(cur_anc_bl);
			for each in ancs:
				if each != node and ancs[each] == cur_anc:
					cur_sis = each;
					nofo[node].append(cur_sis);
					if tree_type == 1:
						cur_sis_bl = getBranchLength(new_tree,cur_sis);
					elif tree_type == 2:
						cur_sis_bl = "NA";
					nofo[node].append(cur_sis_bl);
		else:
			j = 0;
			while j < 4:
				nofo[node].append("");
				j = j + 1;

		tipflag = 1;

		for each in ancs:
			if ancs[each] == node:
				tipflag = 0;
				cur_desc = each;
				nofo[node].append(cur_desc);
				if tree_type == 1:
					cur_desc_bl = getBranchLength(new_tree,cur_desc);
				elif tree_type == 2:
					cur_desc_bl = "NA";
				nofo[node].append(cur_desc_bl);

		if tipflag == 1:
			j = 0;
			while j < 4:
				nofo[node].append("");
				j = j + 1;

		if nofo[node][8] == "":
			nofo[node].append("tip");
		elif nofo[node][0] == None:
			nofo[node].append("root");
		else:
			nofo[node].append("internal");

	##End info retrieval block.

#	for key in nofo:
#		print key + ":" + str(nofo[key]);


	return nofo, new_tree;

#############################################################################
