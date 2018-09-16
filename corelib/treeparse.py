#!/usr/bin/python
#############################################################################
#Functions for retrieving information from newick formatted, rooted phylogenetic trees.
#Gregg Thomas
#Spring 2013-present
#############################################################################

import sys, re

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

def remBranchLength(treestring):
# Removes branch lengths from a tree.

	treestring = re.sub('[)][\d\w<>.eE_:-]+', ')', treestring);
	treestring = re.sub(':[\d.eE-]+', '', treestring);
	#treestring = re.sub('[)][_\d\w<>.eE-]+:[\d.eE-]+', ')', treestring);
	#treestring = re.sub(':[\d.eE-]+', '', treestring);
	#treestring = re.sub('<[\d\w]+>[\d_]+', '', treestring);
	#treestring = re.sub('<[\d\w]+>', '', treestring);
	#treestring = re.sub('[\d_]+', '', treestring);

	return treestring;

#############################################################################

def addBranchLength(tree, treedict):
# Re-writes the branch lengths onto a topology parsed by treeParse.
	for node in treedict:
		if treedict[node][2] == 'root':
			continue;
		if treedict[node][0] != "NA":
			tree = tree.replace(node, node + ":" + treedict[node][0]);
			if treedict[node][3] != "NA":
				tree = tree.replace(node + ":" + treedict[node][0], node + "_" + treedict[node][3] + ":" + treedict[node][0]);
		elif treedict[node][3] != "NA":
			tree = tree.replace(node, node + "_" + treedict[node][3]);
	return tree;

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
		if c_treedict[d][2] != 'tip':
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

def getCladeNode(c_spec, c_treedict):
# This function takes a node in the current tree and the dictionary of the current tree
# (returned by treeparse) and finds all tip labels that are descendants of the current node.
# This is done by getting the direct descendants of the current node with getDesc and then
# recursively calling itself on those descendants.

	clade = [];
	c_desc = getDesc(c_spec, c_treedict);
	for d in c_desc:
		if c_treedict[d][2] != 'tip':
			clade.append(getCladeNode(d, c_treedict));
			# Recursion
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

def pathToRoot(node, tree_dict):
	#ptr = [node];
	ptr = [];
	while tree_dict[node][2] != 'root':
		ptr.append(tree_dict[node][1]);
		node = tree_dict[node][1];
	return ptr;

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
		if treedict[spec][2] == 'root':
			continue;
		curanc = treedict[spec][1];
		ancs[spec].append(curanc);
		while treedict[curanc][2] != 'root':
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

def getSubtree(node, tree):
	subtree = "";
	partree = tree[:tree.index(node)][::-1];
	cp = 0;
	op = 0;
	for c in partree:
		if c == ")":
			cp = cp + 1;
		if c == "(":
			op = op + 1;
		subtree = subtree + c;
		if cp == op:
			break;
	return subtree[::-1];

#############################################################################
def comAnc(spec_list, treedict):
#Given a list of species within the tree and the dictionary returned by treeParse using that tree,
#this function checks whether those species are monophyletic (ie they all share a common ancestor).

	#print spec_list;
	if len(spec_list) > 1:
		if spec_list.count(spec_list[0]) == len(spec_list):
			if treedict[spec_list[0]][2] == 'root':
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
#This function returns a list of the nodes between the current node and the root of the tree.
	if n_treedict[n_spec][2] == 'root':
		return [];

	ancs = []
	curanc = n_treedict[n_spec][1];
	ancs.append(curanc);
	while n_treedict[curanc][2] != 'root':
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
#This function counts the number of internal nodes in a tree.
	num_nodes = 0;
	for node in treedict:
		if treedict[node][2] != 'tip':
			num_nodes = num_nodes + 1;
	return num_nodes;

#############################################################################

def nodeDist(query_node, target_node, treedict):
	dist = float(treedict[query_node][0]);
	ancnode = treedict[query_node][1];
	#print query_node, '-', dist, '-', ancnode, '-', treedict[ancnode][0];
	if ancnode == target_node:
		return dist;
	else:
		dist += nodeDist(ancnode, target_node, treedict);

	return dist;

#############################################################################

def rootedOrNot(treedict):
	num_tips = len([n for n in treedict if treedict[n][2] == 'tip']);
	num_internal = len([n for n in treedict if treedict[n][2] != 'tip']);
	if num_internal != (num_tips - 1):
		return 0;
	elif num_internal == (num_tips - 1):
		return 1;
	else:
		return -1;

#############################################################################

def treeParse(tree, debug=0):
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# New (current) format:
# node:[branch length (if present), ancestral node, node type, node label (if present)]

	tree = tree.strip();
	if tree[-1] != ";":
		tree += ";";
	# Some string handling

	nodes, bl, supports, ancs = {}, {}, {}, {};
	# Initialization of all the tracker dicts

	topology = remBranchLength(tree);

	if debug == 1:
		print "TOPOLOGY:", topology;

	nodes = {};
	for n in topology.replace("(","").replace(")","").replace(";","").split(","):
		nodes[n] = 'tip';
	# nodes = { n : 'tip' for n in topology.replace("(","").replace(")","").replace(";","").split(",") };
	# Retrieval of the tip labels

	if debug == 1:
		print "NODES:", nodes;

	new_tree = "";
	z = 0;
	numnodes = 1;
	while z < (len(tree)-1):
		new_tree += tree[z];
		if tree[z] == ")":
			node_label = "<" + str(numnodes) + ">";
			new_tree += node_label;
			nodes[node_label] = 'internal';
			numnodes += 1;
		z += 1;
	nodes[node_label] = 'root';
	rootnode = node_label;
	# This labels the original tree as new_tree and stores the nodes and their types in the nodes dict

	if debug == 1:
		print "NEW TREE:", new_tree;
		print "TREE:", tree;
		print "NODES:", nodes;
		print "ROOTNODE:", rootnode;
		#sys.exit();
	topo = "";
	z = 0;
	numnodes = 1;
	while z < (len(topology)-1):
		topo += topology[z];
		if topology[z] == ")":
			node_label = "<" + str(numnodes) + ">";
			topo += node_label;
			numnodes += 1;
		z += 1;
	# This labels the topology with the same internal labels

	if debug == 1:
		print "TOPO:", topo;
		print "----------";
		print "TOPOLOGY:", topo;

	for node in nodes:
		if node + node in new_tree:
			new_tree = new_tree.replace(node + node, node);

	# if debug == 1:
	# 	print new_tree;
	# 	sys.exit();

	for node in nodes:
	# One loop through the nodes to retrieve all other info
		if debug == 1:
			print "NODE:", node;

		if nodes[node] == 'tip':
			supports[node] = "NA";
			if node + ":" in tree:
				cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
				cur_bl = cur_bl[0].replace(node + ":", "");
				if debug == 1:
					print "FOUND BL:", cur_bl;
				bl[node] = cur_bl;				
			else:
				bl[node] = "NA";

		elif nodes[node] == 'internal':
			if node + node in new_tree:
				new_tree = new_tree.replace(node + node, node);

			if node + "(" in new_tree or node + "," in new_tree or node + ")" in new_tree:
				if debug == 1:
					print "NO BL OR LABEL";
				supports[node] = "NA";
				bl[node] = "NA";

			elif node + ":" in new_tree:
				supports[node] = "NA";
				cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
				cur_bl = cur_bl[0].replace(node + ":", "");
				if debug == 1:
					print "FOUND BL:", cur_bl;
				bl[node] = cur_bl;								

			else:
				cur_bsl = re.findall(node + "[\d\w<>_*+.Ee-]+:[\d.Ee-]+", new_tree);
				if cur_bsl:
				# If the pattern above is found then the node has both support and branch length
					cur_bs = cur_bsl[0].replace(node, "");
					cur_bs = cur_bs[:cur_bs.index(":")];
					cur_bl = cur_bsl[0].replace(node, "").replace(cur_bs, "").replace(":", "");
					if debug == 1:
						print "FOUND BL AND LABEL:", cur_bl, cur_bs;
					supports[node] = cur_bs;
					bl[node] = cur_bl;
					#new_tree = new_tree.replace(cur_bs, "");
				else:
				# If it is not found then the branch only has a label
					cur_bs = re.findall(node + "[\w*+.<> -]+", new_tree);
					cur_bs = cur_bs[0].replace(node, "");
					if debug == 1:
						print "FOUND LABEL:", cur_bs;
					supports[node] = cur_bs;
					bl[node] = "NA";
					#new_tree = new_tree.replace(cur_bs, "");

		elif nodes[node] == 'root':
			bl[node] = "NA";
			supports[node] = new_tree[new_tree.index(node)+len(node):];
			ancs[node] = "NA";
			continue;

		# Next we get the ancestral nodes. If the node is the root this is set to NA.
		anc_match = re.findall('[(),]' + node, new_tree);

		#if nodes[node] == 'internal':
		#	sys.exit();
		# anc_match = re.findall(node + '[\d:(),]+', new_tree);
		anc_match = re.findall(node, topo);
		if debug == 1:
			print "ANC MATCH:", anc_match;

		anc_tree = new_tree[new_tree.index(anc_match[0]):][1:];
		# Ancestral labels are always to the right of the node label in the text of the tree, so we start our scan from the node label

		if debug == 1:
			print "NODE:", node;
			print "ANC_MATCH:", anc_match;
			print "ANC_TREE:", anc_tree;
			
		cpar_count = 0;
		cpar_need = 1;

		for i in range(len(anc_tree)):
		# We find the ancestral label by finding the ) which matches the nesting of the number of ('s found
			if anc_tree[i] == "(":
				cpar_need = cpar_need + 1;
			if anc_tree[i] == ")" and cpar_need != cpar_count:
				cpar_count = cpar_count + 1;
			if anc_tree[i] == ")" and cpar_need == cpar_count:
				anc_tree = anc_tree[i+1:];
				ancs[node] = anc_tree[:anc_tree.index(">")+1];
				break;

		if debug == 1:
			print "FOUND ANC:", ancs[node];
			print "---";
	nofo = {};
	for node in nodes:
		nofo[node] = [bl[node], ancs[node], nodes[node], supports[node]];
	# Now we just restructure everything to the old format for legacy support

	if debug == 1:
	# Debugging options to print things out
		print("\ntree:\n" + tree + "\n");
		print("new_tree:\n" + new_tree + "\n");
		print("topology:\n" + topo + "\n");
		print("nodes:");
		print(nodes);
		print
		print("bl:");
		print(bl);
		print
		print("supports:");
		print(supports);
		print
		print("ancs:");
		print(ancs);
		print
		print("-----------------------------------");
		print
		print("nofo:");
		print(nofo);
		print

	return nofo, topo, rootnode;

#############################################################################















## These are BOTH now old tree parsing functions. Just keeping them around until I manage to replace everything...

#############################################################################
def treeParseNew(tree, tree_type):
#The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
#dictionary with usable info about the tree in the following format:
#New (current) format:
#node:[branch length, ancestral node, ancestral branch length, node type, node labels (if present)]
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
	supports = {};

	while z < (len(tree)-1):
		if tree_type == 1:
			if tree[z] == ":" and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
			elif tree[z] == ":":
				tmp_str = tree[:z];
				if tmp_str.rfind(")") > tmp_str.rfind(","):
					new_node = "<" + str(numnodes) + ">";

					supports[new_node] = tmp_str[tmp_str.rfind(")")+1:];

					new_tree = new_tree[:new_tree.rfind(")")+1] + new_node;
					numnodes = numnodes + 1;

		if tree_type == 2:
			if (tree[z] == "," or tree[z] == ")") and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		new_tree = new_tree + tree[z];
		z = z + 1;

	if new_tree[-1] not in [")",">"]:
		if new_tree.rfind(")") > new_tree.rfind(">"):
			last_char = ")";
		else:
			last_char = ">";

		new_tree = new_tree[:new_tree.rfind(last_char)+1];

	if new_tree[-1] == ")":
		rootnode = "<" + str(numnodes) + ">"
		new_tree = new_tree + rootnode;
	else:
		rootnode = new_tree[new_tree.rfind(")")+1:];

	##This first block labels all internal nodes with the format <#>
	
	# print tree;
	# print new_tree;
	# print supports;
	# print "-----------------------------------";

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
						#if a == (len(new_tree)-5):
						#	curanc = new_tree[a+1:];
						if new_tree[a+1:].find(":") == -1:
							#curanc = new_tree[len(new_tree)-5:];
							curanc = new_tree[new_tree.rfind(")")+1:]
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
	#for key in ancs:
	#	print key + ":", ancs[key]
	#print "---------";
	#sys.exit()

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

		if nofo[node][2] != 'tip' and supports != {} and node in supports:
			nofo[node].append(supports[node]);

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
def treeParseOld(tree, tree_type):
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