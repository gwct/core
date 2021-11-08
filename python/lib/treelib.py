#############################################################################
# Tree modules -- used by tree.py
# Gregg Thomas
# August 2017
#############################################################################

import core, sys, os, subprocess, treeparse as tp, re
from collections import defaultdict

#############################################################################

def treeSep(infile, out_prefix, outdir):
# This function takes an input file and separates all trees in it to their own files.
	if not out_prefix:
		out_prefix = os.path.splitext(os.path.basename(infile))[0];
	# If no prefix is specified, one is chosen based on the input file name.

	pathfile = open(os.path.join(outdir, "tree-sep-paths.txt"), "w");
	# Creates a file with a list of paths to each gene tree file. Thought I would use this for Notung, but maybe not.

	num_lines, num_trees, line_skip = 0,0,[];
	for line in open(infile):
		num_lines += 1;
		tree = line.strip();
		try:
			td, t, r = tp.treeParse(line);
		except:
			line_skip.append(str(num_lines-1));
			continue;
		num_trees += 1;
		# Check to make sure each line can be read as a Newick string.

		outfilename = os.path.join(outdir, out_prefix + "-" + str(num_lines) + ".tre");
		pathfile.write(os.path.abspath(outfilename) + "\n");
		with open(outfilename, "w") as treefile:
			treefile.write(tree);
		# Write the valid Newick string to its own output file.
	pathfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print(str(num_lines) + " lines in file.");
	if line_skip != []:
		print("The following " + str(len(line_skip)) + " line(s) were skipped because they couldn't be read as Newick formatted trees: " + ",".join(line_skip));
	print(str(num_trees) + " trees separated!");
	print("=======================================================================");

#############################################################################

def treeJoin(infiles, outfilename):
# This function takes an input directory and reads all files. It puts any Newick strings
# it finds into the output file.
	num_files, num_read, num_trees, tre_skip, parse_skip = 0,0,0,[],[];
	with open(outfilename, "w") as treefile:
		for infile in infiles:
			num_files += 1;
			for line in open(infile):
				line = line.strip();
				try:
					td, t, r = tp.treeParse(line);
				except:
					if infile not in parse_skip:
						parse_skip.append(infile);
					continue;
				num_trees += 1;
				treefile.write(line + "\n");
				# Check if each line in the current file is a Newick string and, if so, write it to 
				# the output file.

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print(str(num_files) + " total files.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " file(s) were skipped because they couldn't be read as tree files: " + ",".join([os.path.basename(f) for f in tre_skip]));
	print(str(num_trees) + " trees joined.");
	if parse_skip != []:
		print("The following " + str(len(parse_skip)) + " file(s) had lines that couldn't be read as trees and were skipped: " + ",".join([os.path.basename(f) for f in parse_skip]));
	print("=======================================================================");

#############################################################################

def labelTree(infiles, tree_flag, outfilename):
# This function takes either a tree string or a file with multiple trees and
# simply puts integer labels on the internal nodes in the format <#>
	if tree_flag:
		td, tree, r = infiles;
		labeled_tree = tp.addBranchLength(tree, td);
		print("\n" + labeled_tree + "\n");
		sys.exit();
	# For a tree string, just print the labeled tree to the screen

	infile = infiles[0];
	num_lines, num_trees, line_skip = 0,0,[];
	with open(outfilename, "w") as treefile:
		for line in open(infile):
			num_lines += 1;
			try:
				td, tree, r = tp.treeParse(line);
			except:
				line_skip.append(str(num_lines));
				continue;
			num_trees += 1;
			# for each line in the file, check to make sure it is a Newick string.

			labeled_tree = tp.addBranchLength(tree, td);
			treefile.write(labeled_tree + "\n");
			# Label the tree and write to the output file.

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print(str(num_lines) + " total lines.");
	print(str(num_trees) + " trees labeled.");
	if line_skip != []:
		print("The following " + str(len(line_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join([os.path.basename(f) for f in line_skip]));
	print("=======================================================================");	

#############################################################################

def rootCheck(infiles, tree_flag, outfilename):
# This function takes an input Newick string or file and checks if the trees are rooted or not.
	if tree_flag:
		td, tree, r = infiles;
		rooted = tp.rootedOrNot(td);
		if rooted == 0:
			print("\nUnrooted!\n")
		elif rooted == 1:
			print("\nRooted!\n");
		else:
			print("\n???\n");
		sys.exit();
	# If the input is a Newick string, simply print the result to the screen.

	infile = infiles[0];
	num_lines, num_trees, line_skip, num_unroot, num_rooted, num_undet = 0,0,[],0,0,0;
	with open(outfilename, "w") as treefile:
		for line in open(infile):
			num_lines += 1;
			try:
				td, tree, r = tp.treeParse(line);
			except:
				line_skip.append(str(num_lines));
				treefile.write("**Skipped\n");
				continue;
			num_trees += 1;
			# Check to make sure each line in the file can be read as a Newick string. If not, skip it.

			rooted = tp.rootedOrNot(td);
			# This function works on the basis that an unrooted tree will be read as a trifurcation at the root by my
			# tree parser. If that is the case, it will have one fewer internal nodes than expected. 
			# This means that this will only work on trees that have only bifurcating internal nodes.

			if rooted == 0:
				num_unroot += 1;
				treefile.write("Unrooted!\n")
			if rooted == 1:
				num_rooted += 1;
				treefile.write("Rooted!\n");
			else:
				num_undet += 1;
				treefile.write("???\n");

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print(str(num_lines) + " total lines.");
	if line_skip != []:
		print("The following " + str(len(line_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(line_skip));
	print(str(num_trees) + " trees checked.");
	print(str(num_rooted) + " rooted trees.");
	print(str(num_unroot) + " unrooted trees.");
	print(str(num_undet) + " undetermined trees.");
	print("=======================================================================");	

#############################################################################

def rootTrees(infiles, tree_flag, outgroup, outfilename):
# This function relies on Newick Utilities (NU) to root trees at a specified outgroup.
	import re
	tmpfilename = "tmp9825xyz-t-m-p.tmp";
	# NU can only read trees from a file, so I make a temporary one.
	try:
		outgroup = outgroup.split(",");
	except:
		sys.exit(core.errorOut(26, "-outgroup entered incorrectly! Should be comma delimited list of tip labels."));
	# Check to make sure the outgroups were entered correctly.

	if tree_flag:
		td, tree, r, tree_string = infiles;
		lca, monophyletic = tp.LCA(outgroup, td);
		if monophyletic == 0:
			sys.exit(core.errorOut(27, "Your outgroup labels (-outgroup) must be monophyletic!"));
		# Specified outgroups must be monophyletic.
		
		with open(tmpfilename, "w") as tmpfile:
			tmpfile.write(tree_string);
		nw_cmd = "nw_reroot " + tmpfilename + " " + " ".join(outgroup);
		print("\n----Rooted tree----");
		os.system(nw_cmd);
		print()
		# The NU call with nw_reroot.
		sys.exit();
	# If the input is a Newick string, just print the output to the screen.

	tmpfilename_2 = "tmp25xzgz-t-m-p.tmp"
	# To retrieve output from NU I use another temporary file.
	infile = infiles[0];
	num_lines, num_trees, non_mono, line_skip = 0,0,[],[];
	with open(outfilename, "w") as treefile:
		for line in open(infile):
			num_lines += 1;
			tid = False;
			if len(line.split("\t")) > 1:
				tid, line = line.strip().split("\t");
			try:
				td, tree, r = tp.treeParse(line);
			except:
				line_skip.append(str(num_lines));
				treefile.write("**Skipped - couldn't read as Newick string\n");
				continue;
			if not all(s in td for s in outgroup):
				line_skip.append(str(num_lines));
				treefile.write("**Skipped - not all outgroups in tree.\n");
				continue;

			lca, monophyletic = tp.LCA(outgroup, td);
			if monophyletic == 0:
				non_mono.append(str(num_lines));
				continue;
			num_trees += 1;
			# Check to make sure each line is a Newick string and that the outgroups are monophyletic in that tree.

			with open(tmpfilename, "w") as tmpfile:
				tmpfile.write(line);
			nw_cmd = "nw_reroot " + tmpfilename + " " + " ".join(outgroup) + " > " + tmpfilename_2;
			os.system(nw_cmd);
			# The NU call with nw_reroot.

			rooted_tree = open(tmpfilename_2, "r").read().strip();
			if tid:
				rooted_tree = tid + "\t" + rooted_tree;
			treefile.write(rooted_tree + "\n");
			# Getting the output from the tmp file.

	os.system("rm " + tmpfilename);
	os.system("rm " + tmpfilename_2);
	# Remove the tmp files.

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print(str(num_lines) + " total lines.");
	if line_skip != []:
		print("The following " + str(len(line_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(line_skip));
	if non_mono != []:
		print("The following " + str(len(non_mono)) + " lines did not contain monophyletic outgroups: " + ",".join(non_mono));
	print(str(num_trees) + " trees rooted.");
	print("=======================================================================");

#############################################################################

def rootTreesBest(infiles, tree_flag, outgroup, outfilename):
# This function relies on Newick Utilities (NU) to root trees at a specified outgroup.
	import re
	tmpfilename = "tmp9825xyz-t-m-p.tmp";
	# NU can only read trees from a file, so I make a temporary one.
	try:
		outgroup = outgroup.split(",");
	except:
		sys.exit(core.errorOut(26, "-outgroup entered incorrectly! Should be comma delimited list of tip labels."));
	# Check to make sure the outgroups were entered correctly.

	tmpfilename_2 = "tmp25xzgz-t-m-p.tmp"
	# To retrieve output from NU I use another temporary file.
	infile = infiles[0];
	num_lines, num_trees, non_mono, line_skip = 0,0,[],[];
	with open(outfilename, "w") as treefile:
		for line in open(infile):
			alnfilename, line = line.strip().split("\t");
			num_lines += 1;
			try:
				td, tree, r = tp.treeParse(line);
			except:
				line_skip.append(str(num_lines));
				treefile.write("**Skipped - couldn't read as Newick string\n");
				continue;
			if not all(s in td for s in outgroup):
				line_skip.append(str(num_lines));
				treefile.write("**Skipped - not all outgroups in tree.\n");
				continue;

			lca, monophyletic = tp.LCA(outgroup, td);
			if monophyletic == 0:
				non_mono.append(str(num_lines));
				continue;
			num_trees += 1;
			# Check to make sure each line is a Newick string and that the outgroups are monophyletic in that tree.

			with open(tmpfilename, "w") as tmpfile:
				tmpfile.write(line);
			nw_cmd = "nw_reroot " + tmpfilename + " " + " ".join(outgroup) + " > " + tmpfilename_2;
			os.system(nw_cmd);
			# The NU call with nw_reroot.

			rooted_tree = open(tmpfilename_2, "r").read().strip();
			treefile.write(alnfilename + "\t" + rooted_tree + "\n");
			# Getting the output from the tmp file.

	os.system("rm " + tmpfilename);
	os.system("rm " + tmpfilename_2);
	# Remove the tmp files.

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print(str(num_lines) + " total lines.");
	if line_skip != []:
		print("The following " + str(len(line_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(line_skip));
	if non_mono != []:
		print("The following " + str(len(non_mono)) + " lines did not contain monophyletic outgroups: " + ",".join(non_mono));
	print(str(num_trees) + " trees rooted.");
	print("=======================================================================");

#############################################################################

def flightOfTheConcordance(infiles, tree_flag, genefilename, count_tops):
# This function calculates concordance factors for each node in a species tree given a
# set of singly-copy gene trees.
	if tree_flag:
		sinfo, stree, sroot, sfilename = infiles;
	else:
		try:
			sinfo, stree, sroot = tp.treeParse(open(infiles[0], "r").read());
		except:
			sys.exit(core.errorOut(28, "Could not read species tree (-s) as a Newick tree!"));
		# If the input species tree was a file check to make sure it contains a valid Newick tree.	

	stips = [node for node in sinfo if sinfo[node][1] == 'tip'];
	sclades = { node : set(tp.getClade(node, sinfo)) for node in sinfo if sinfo[node][2] != 'tip' };
	# Get the tips and clades (tip nodes) for each internal node in the species tree.

	#node_counts = defaultdict(float);
	node_counts = { node : 0.0 for node in sinfo if sinfo[node][2] != 'tip' };
	num_lines, tre_skip, sc_skip = 0, [], [];
	total_trees = 0.0;
	if count_tops:
		tops, top_counts, top_trees = [], [], [];

	for line in open(genefilename):
		num_lines += 1;
		try:
			ginfo, gtree, groot = tp.treeParse(line);
		except:
			tre_skip.append(str(num_lines));
			continue;
		# Check if each line in the genetrees file is a Newick string.

		if count_tops:
			gclade = set([frozenset(tp.getClade(node, ginfo)) for node in ginfo if ginfo[node][2] != 'tip']);
			if gclade in tops:
				topind = tops.index(gclade);
				top_counts[topind] += 1;
			else:
				tops.append(gclade);
				top_counts.append(1);
				top_trees.append(gtree);

		gtips = [node for node in ginfo if ginfo[node][1] == 'tip'];
		if set(gtips) != set(stips) or len(gtips) != len(stips):
			sc_skip.append(str(num_lines));
			continue;
		# Check to make sure the tips are identical in the current gene tree and the species tree.

		gclades = [ set(tp.getClade(node, ginfo)) for node in ginfo if ginfo[node][2] != 'tip' ];
		# Get the clades (tip nodes) for each internal node in the current gene tree.

		for node in sclades:
			if sclades[node] in gclades:
				node_counts[node] += 1.0;
			if node == sroot and sclades[node] not in gclades:
				print(len(ginfo));
				print(gtree);
		# Check if each species tree clade is in the gene tree.
		total_trees += 1.0;


	print("\n" + core.getTime() + " Done!");

	if count_tops:
		print("\n----Topology counts----");
		tops_dict = {};
		for x in range(len(tops)):
			tops_dict[top_trees[x]] = top_counts[x];
		for item in sorted(list(tops_dict.items()), key=lambda x: x[1], reverse=True):
			print(item[0], item[1]);
		print(len(tops_dict), "total topologies found");

	print("\n----Concordance factor nodes----");
	stree = tp.addBranchLength(stree, sinfo);
	for node in node_counts:
		cf = round(node_counts[node]/total_trees,2)
		print(node, cf);
		if sinfo[node][3] not in  ["NA",'']:
			stree = stree.replace(node, node + "_" + str(cf));
		else:
			stree = stree.replace(node, node+ "_" + str(cf));

	if stree[-1] == "_":
		stree = stree[:-1];

	print("\n----Concordance factor tree----");
	print(stree);
	print()

	print("-----");
	print(str(num_lines) + " total lines in gene tree file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	if sc_skip != []:
		print("The following " + str(len(sc_skip)) + " lines were skipped because they didn't have the same number of nodes as the species tree: " + ",".join(sc_skip));
	print(str(total_trees) + " trees read.");
	print("=======================================================================");

#############################################################################

def countTips(infile):
# This function counts all unique tip labels given a set of trees.

	num_lines, num_trees, tre_skip = 0, 0, [];

	tips = defaultdict(int);

	for line in open(infile):
		num_lines += 1;
		try:
			td, tree, root = tp.treeParse(line);
		except:
			tre_skip.append(str(num_lines));
			continue;
		num_trees += 1;
		# Check if each line in the genetrees file is a Newick string.

		for node in td:
			if td[node][2] == 'tip':
				tips[node] += 1;
		# Iterate the dictionary for each tip in the current tree.

	maxlen = 0;
	for tip in tips:
		if len(tip) > maxlen:
			maxlen = len(tip);
	# Get the length of the tip labels for some nicer printing.

	print("\n" + core.getTime() + " Done!");
	print("\n----Tip counts----");
	pad = maxlen + 2;
	for tip in tips:
		print(core.spacedOut(tip, pad), tips[tip]);
	# Print the tip labels and counts.
	print("\n-----");
	print(str(num_lines) + " total lines in input file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	print(str(num_trees) + " trees read.");
	print("=======================================================================");

#############################################################################

def countClade(infile, clade):
# This function counts all unique tip labels given a set of trees.

	num_lines, num_trees, tre_skip, clade_count = 0, 0, [], 0;

	for line in open(infile):
		num_lines += 1;
		try:
			td, tree, root = tp.treeParse(line);
		except:
			tre_skip.append(str(num_lines));
			continue;
		num_trees += 1;
		# Check if each line in the genetrees file is a Newick string.

		for node in td:
			if td[node][2] != 'tip':
				cur_clade = set(tp.getClade(node, td));
				if cur_clade == clade:
					clade_count += 1;
		# Iterate the dictionary for the clade in the current tree.

	print("\n" + core.getTime() + " Done!");
	print("\n----Clade counts----");
	print("# of trees with clade:\t", clade_count);
	# Print the # of trees containing the clade.
	print("\n-----");
	print(str(num_lines) + " total lines in input file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	print(str(num_trees) + " trees read.");
	print("=======================================================================");

#############################################################################
def relabelTips(infile, labels, mode, delim, output):
# This function takes a file with many trees and searches the tip labels to match
# strings for replacements.
	try:
		labels = {l.split(",")[0] : l.split(",")[1] for l in labels.split(" ")};
	except:
		sys.exit(core.errorOut(29, "-labels was not input correctly! Format 'oldlabel1,newlabel1 oldlabel2,newlabel2'"));
	# Check to make sure the labels were input properly by the user.

	if delim == 'space':
		delim = ' ';

	pad = 20;
	print("\n---------Relabeling tips---------");
	print(core.spacedOut("Old label contains", pad), "| New label");
	print("---------------------------------");
	for old in labels:
		print(core.spacedOut(old, pad), "| " + labels[old]);
	# Some nice printing of the labels.

	num_lines, num_trees, tre_skip = 0, 0, [];
	with open(output, "w") as outfile:
		for line in open(infile):
			line = line.strip();
			num_lines += 1;
			try:
				td, tree, root = tp.treeParse(line);
			except:
				tre_skip.append(str(num_lines));
				continue;
			num_trees += 1;
			# Check if each line in the genetrees file is a Newick string.

			for node in td:
				for old in labels:
					if old in node:
						if mode == 1:
							line = line.replace(node, labels[old]);
						if mode == 2:
							line = line.replace(node, labels[old] + delim + node);
						if mode == 3:
							line = line.replace(node, node + delim + labels[old]);
			# Check and relabel every tip in the current tree.
			# mode == 1 : replace old label with new label
			# mode == 2 : put new label on beginning of old label
			# mode == 3 : put new label at end of old label

			outfile.write(line + "\n");
			# For each node in the tree, check if it contains the text of one of the labels to replace.
			# If so, replace it.

	print("---------------------------------");
	print("\n" + core.getTime() + " Done!");
	print(str(num_lines) + " total lines in input file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	print(str(num_trees) + " trees read.");
	print("=======================================================================");

#############################################################################

def rmLabel(infile, mode, outfilename, best_flag=False):
# Takes a file with many trees and removes internal node labels and/or branch lengths (depending on mode).
	num_lines, num_trees, tre_skip, parse_skip = 0,0,[],[];
	with open(outfilename, "w") as treefile:
		for line in open(infile):
			num_lines +=1;
			line = line.strip();
			if best_flag:
				title, line = line.split("\t");
			try:
				td, out_tree, r = tp.treeParse(line);
			except:
				if infile not in parse_skip:
					parse_skip.append(infile);
				continue;
			num_trees += 1;
			# Check if each line in the genetrees file is a Newick string.

			if mode == 1:
				out_tree = re.sub('\)[\d\w]+:', '):', line);
			if mode == 2:
				out_tree = re.sub(':[\d.eE-]+', '', line);
			if mode == 3:
				out_tree = re.sub('<[\d]+>', '', out_tree);
			# mode == 1 : remove internal node labels only
			# mode == 2 : remove branch lengths only
			# mode == 3 : remove both internal node labels and branch lengths

			if best_flag:
				treefile.write(title + "\t" + out_tree + "\n");
			else:
				treefile.write(out_tree + "\n");
			# Write the edited tre to the output file.

	print("\n-----");
	print("\n" + core.getTime() + " Done!");
	print(str(num_lines) + " total lines in input file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	print(str(num_trees) + " trees read.");
	print("=======================================================================");

#############################################################################

def scaleBL(infile, op, factor, outfilename):
# Takes a file with many trees and removes internal node labels and/or branch lengths (depending on mode).
	num_lines, num_trees, tre_skip, parse_skip = 0,0,[],[];
	with open(outfilename, "w") as treefile:
		for line in open(infile):
			num_lines +=1;
			line = line.strip();
			try:
				td, out_tree, r = tp.treeParse(line);
			except:
				if infile not in parse_skip:
					parse_skip.append(infile);
				continue;
			num_trees += 1;
			# Check if each line in the genetrees file is a Newick string.

			for node in td:
				if node == r:
					continue;
				old_bl = float(td[node][0]);
				if op == "/":
					new_bl = old_bl / factor;
				elif op == "*":
					new_bl = old_bl * factor;
				elif op == "+":
					new_bl = old_bl + factor;
				elif op == "-":
					new_bl = old_bl - factor;
				td[node][0] = str(new_bl);

			out_tree = tp.addBranchLength(out_tree, td);
			treefile.write(out_tree + "\n");
			# Write the edited tre to the output file.

	print("\n-----");
	print("\n" + core.getTime() + " Done!");
	print(str(num_lines) + " total lines in input file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	print(str(num_trees) + " trees read.");
	print("=======================================================================");

#############################################################################

def robF(infiles, tree_flag, genefilename, raxpath, outfile):
# This function calls RAxML to calculate Robinson-Foulds distances for each gene tree to the species tree.
	import re
	raxlog = os.path.splitext(outfile)[0] + ".raxlog";

	if tree_flag:
		sinfo, stree, sroot, sfilename = infiles;
	else:
		try:
			sinfo, stree, sroot = tp.treeParse(open(infiles[0], "r").read());
		except:
			sys.exit(core.errorOut(30, "Could not read species tree (-s) as a Newick tree!"));
		# If the input species tree was a file check to make sure it contains a valid Newick tree.	

	stips = [node for node in sinfo if sinfo[node][1] == 'tip'];
	# Get the tips in the species tree.

	stree = re.sub("<[\d]+>", "", stree);
	# Remove node labels from species tree.

	num_lines, tre_skip, sc_skip, raxfail = 0, [], [], [];
	rfs, wrfs = [], [];
	total_trees = 0.0;
	tmpfile = "tmp74ghr2.tmp";
	rfoutfile = "RAxML_RF-Distances.RFtmp7f";

	with open(outfile, "w") as out:
		for line in open(genefilename):
			num_lines += 1;
			try:
				ginfo, gtree, groot = tp.treeParse(line);
			except:
				out.write("Couldn't read as Newick string -- Skipping.\n");
				tre_skip.append(str(num_lines));
				continue;
			# Check if each line in the genetrees file is a Newick string.

			gtips = [node for node in ginfo if ginfo[node][1] == 'tip'];
			if set(gtips) != set(stips) or len(gtips) != len(stips):
				out.write("Tip labels not identical to species tree -- Skipping.\n");
				sc_skip.append(str(num_lines));
				continue;
			# Check to make sure the tips are identical in the current gene tree and the species tree.

			gtree = re.sub("<[\d]+>", "", gtree);
			# Remove node labels from gene tree.

			with open(tmpfile, "w") as tfile:
				tfile.write(stree + ";\n");
				tfile.write(gtree + ";");

			raxcall = raxpath + " -m PROTGAMMAJTTF -z " + tmpfile + " -f r -n RFtmp7f >> " + raxlog;
			#print num_lines, raxcall;
			os.system(raxcall);

			if os.path.exists(rfoutfile):
				curout = open(rfoutfile, "r").read().strip().split(" ");
				cur_rf, cur_wrf = curout[2], curout[3];
				rfs.append(float(cur_rf));
				wrfs.append(float(cur_wrf));
				out.write(cur_rf + " " + cur_wrf + "\n");
			else:
				out.write("RAxML failed -- Skipping.\n");
				raxfail.append(str(num_lines));

			os.system("rm " + rfoutfile);
			os.system("rm RAxML_info.RFtmp7f");
			total_trees += 1.0;

	print("-----");
	print("Average RF distance for all gene trees:", float(sum(rfs)) / float(len(rfs)));
	print("Average weighted RF distance for all gene trees:", round(float(sum(wrfs)) / float(len(wrfs)),3));
	os.system("rm " + tmpfile);
	print("-----");
	print(str(num_lines) + " total lines in gene tree file.");
	if tre_skip != []:
		print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
	if sc_skip != []:
		print("The following " + str(len(sc_skip)) + " lines were skipped because they didn't have the same number of nodes as the species tree: " + ",".join(sc_skip));
	if sc_skip != []:
		print("The following " + str(len(sc_skip)) + " lines were skipped because they didn't have the same number of nodes as the species tree: " + ",".join(sc_skip));
	print(str(total_trees) + " trees successfully read and calculated RF distances.");
	print("=======================================================================");