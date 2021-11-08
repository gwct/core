#!/usr/bin/python
#############################################################################
# Functions for retrieving information from newick formatted, phylogenetic trees.
# Gregg Thomas
# Spring 2013-present
#############################################################################

import sys, re

#############################################################################
def getBranchLength(bltree, spec_label):
# Returns the branch length of a species given a newick formatted tree. Used by treeParse.

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

				for a in range(len(indcheck)):
					if indcheck[a] == -1:
						indcheck[a] = 10000;

				curbranch = bltree[d+1:min(indcheck)];
				return curbranch;
		d = d + 1;
	startind = d;

#############################################################################

def remBranchLength(treestring):
# Removes branch lengths from a tree.
	treestring = re.sub('[)][\d\w<>/.eE_:-]+', ')', treestring);
	treestring = re.sub(':[\d.eE-]+', '', treestring);
	return treestring;

#############################################################################

def remNodeLabel(treestring):
# Removes node labels from a tree.
	treestring = re.sub('[)][\d\w<>.eE_-]+', ')', treestring);
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
# Gets all nodes between provided node and root of the tree
	ptr = [];
	while tree_dict[node][2] != 'root':
		ptr.append(tree_dict[node][1]);
		node = tree_dict[node][1];
	return ptr;

#############################################################################

def LCA(spec_list, treedict):
# Given a list of nodes, this function finds the least common ancestor of them,
# and tells whether the nodes provided form a monophyletic clade.

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

	intersect_anc = set.intersection(*list(map(set, list(ancs.values()))))
	lcp = [t for t in list(ancs.values())[0] if t in intersect_anc]

	#lcp = sorted(set.intersection(*map(set, ancs.values())), key=lambda x: ancs.values()[0].index(x))
	monophyletic = False;
	if set(getClade(lcp[0],treedict)) == set(spec_list):
		monophyletic = True;

	return lcp[0], monophyletic;

#############################################################################

def getSubtree(node, tree):
# Gets the subtree string at a given node from a labeled tree string from treeParse

    subtree = "";
    # Initialize the subtree string

    partree = tree[:tree.index(node)][::-1]
    # Slice the tree string at the index of the node label and reverse it

    cp = 0;
    op = 0;
    # Counts of closing an opening parentheses. The subtree will be complete when they
    # are equal

    for c in partree:
        if c == ")":
            cp = cp + 1;
        if c == "(":
            op = op + 1;
        subtree = subtree + c;
        if cp == op:
            break;
    # Loop through every character in the sliced and reversed tree, counting parentheses and adding
    # charcters to the subtree string one at a time. Stop when the number of closing parentheses equals
    # the number of opening

    return subtree[::-1];
    # Return the reverse of the subtree, which is the original orientation

#############################################################################
def comAnc(spec_list, treedict):
# Given a list of species within the tree and the dictionary returned by treeParse using that tree,
# this function checks whether those species are monophyletic (ie they all share a common ancestor).
# Similar to LCA

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
		if list(ancdict.values()).count(ancdict[b]) > 1 and ancdict[b] not in new_list:
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
# This function returns a list of the nodes between the current node and the root of the tree.

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
# Relabels species to match the labels in the tree. (i5k)
	for node in t_d:
		if s in node:
			s = node;
	return s;

#############################################################################

def numTips(treedict):
# This function counts the number of tips in a tree.
	num_nodes = 0;
	for node in treedict:
		if treedict[node][2] == 'tip':
			num_nodes = num_nodes + 1;
	return num_nodes;

#############################################################################

def numInternal(treedict):
# This function counts the number of internal nodes in a tree.
	num_nodes = 0;
	for node in treedict:
		if treedict[node][2] != 'tip':
			num_nodes = num_nodes + 1;
	return num_nodes;

#############################################################################

def nodeDist(query_node, target_node, treedict):
# This gets the distance between two nodes in a tree in terms of their branch length.

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
# Tells whether a tree is rooted or not based on the number of nodes in the
# tree.

	num_tips = len([n for n in treedict if treedict[n][2] == 'tip']);
	num_internal = len([n for n in treedict if treedict[n][2] != 'tip']);
	if num_internal != (num_tips - 1):
		return False;
	elif num_internal == (num_tips - 1):
		return True;
	else:
		return -1;

#############################################################################

def ultrametricOrNot(treedict, root):
# Tells whether a tree is ultra-metric or not based on calls to nodeDist().

	root_dists = [];
	tips = [ n for n in treedict if treedict[n][2] == 'tip' ];
	root_dists = [ round(nodeDist(node, root, treedict), 3) for node in tips ];
	#print(root_dists);
	if root_dists.count(root_dists[0]) == len(root_dists):
		return True;
	else:
		return False;		

#############################################################################

def treeParse(tree, debug=False):
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
		print("TOPOLOGY:", topology);

	##########

	nodes = {};
	for n in topology.replace("(","").replace(")","").replace(";","").split(","):
		nodes[n] = 'tip';
	# Retrieval of the tip labels

	if debug == 1:
		print("NODES:", nodes);

	##########

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

	if debug:
		print("NEW TREE:", new_tree);
		print("TREE:", tree);
		print("NODES:", nodes);
		print("ROOTNODE:", rootnode);
	
	##########

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

	if debug:
		print("TOPO:", topo);

	##########

	for node in nodes:
		if node + node in new_tree:
			new_tree = new_tree.replace(node + node, node);

	##########

	for node in nodes:
	# One loop through the nodes to retrieve all other info
		if debug:
			print("NODE:", node);

		if nodes[node] == 'tip':
			supports[node] = "NA";
			if node + ":" in tree:
				cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
				cur_bl = cur_bl[0].replace(node + ":", "");
				if debug:
					print("FOUND BL:", cur_bl);
				bl[node] = cur_bl;				
			else:
				bl[node] = "NA";

		elif nodes[node] == 'internal':
			if node + node in new_tree:
				new_tree = new_tree.replace(node + node, node);

			if node + "(" in new_tree or node + "," in new_tree or node + ")" in new_tree:
				if debug:
					print("NO BL OR LABEL");
				supports[node] = "NA";
				bl[node] = "NA";

			elif node + ":" in new_tree:
				supports[node] = "NA";
				cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
				cur_bl = cur_bl[0].replace(node + ":", "");
				if debug:
					print("FOUND BL:", cur_bl);
				bl[node] = cur_bl;								

			else:
				cur_bsl = re.findall(node + "[\d\w<>_*+.Ee/-]+:[\d.Ee-]+", new_tree);
				if cur_bsl:
				# If the pattern above is found then the node has both support and branch length
					cur_bs = cur_bsl[0].replace(node, "");
					cur_bs = cur_bs[:cur_bs.index(":")];
					cur_bl = cur_bsl[0].replace(node, "").replace(cur_bs, "").replace(":", "");
					if debug:
						print("FOUND BL AND LABEL:", cur_bl, cur_bs);
					supports[node] = cur_bs;
					bl[node] = cur_bl;
					#new_tree = new_tree.replace(cur_bs, "");
				else:
				# If it is not found then the branch only has a label
					cur_bs = re.findall(node + "[\w*+.<> -]+", new_tree);
					cur_bs = cur_bs[0].replace(node, "");
					if debug:
						print("FOUND LABEL:", cur_bs);
					supports[node] = cur_bs;
					bl[node] = "NA";
					#new_tree = new_tree.replace(cur_bs, "");

		elif nodes[node] == 'root':
			bl[node] = "NA";
			supports[node] = new_tree[new_tree.index(node)+len(node):];
			ancs[node] = "NA";
			continue;
		# Next we get the ancestral nodes. If the node is the root this is set to NA.

		##########

		anc_match = re.findall('[(),]' + node, new_tree);

		anc_match = re.findall(node, topo);
		if debug:
			print("ANC MATCH:", anc_match);

		##########

		anc_tree = new_tree[new_tree.index(anc_match[0]):][1:];
		# Ancestral labels are always to the right of the node label in the text of the tree, so we start our scan from the node label

		if debug:
			print("NODE:", node);
			print("ANC_MATCH:", anc_match);
			print("ANC_TREE:", anc_tree);

		##########
			
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

		if debug:
			print("FOUND ANC:", ancs[node]);
			print("---");

		##########

	####################

	nofo = {};
	for node in nodes:
		nofo[node] = [bl[node], ancs[node], nodes[node], supports[node]];
	# Now we just restructure everything to the old format for legacy support

	if debug:
	# Debugging options to print things out
		print(("\ntree:\n" + tree + "\n"));
		print(("new_tree:\n" + new_tree + "\n"));
		print(("topology:\n" + topo + "\n"));
		print("nodes:");
		print(nodes);
		print()
		print("bl:");
		print(bl);
		print()
		print("supports:");
		print(supports);
		print()
		print("ancs:");
		print(ancs);
		print()
		print("-----------------------------------");
		print()
		print("nofo:");
		print(nofo);
		print()

	return nofo, topo, rootnode;

#############################################################################

