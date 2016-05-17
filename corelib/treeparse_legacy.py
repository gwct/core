import sys, re
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

		if nofo[node][3] != 'tip' and supports != {} and node in supports:
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