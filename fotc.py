#!/usr/bin/python
#############################################################################
# Given a set of single-copy gene trees and a corresponding species tree this
# script calculates concordance factors for all nodes in the species tree.
#
# Dependencies: core
#
# Gregg Thomas, Summer 2017
#############################################################################

import sys, os
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core, treeparse as tr
from collections import defaultdict

if len(sys.argv) != 3 or "-h" in sys.argv:
	print "\nUsage:\t$ python fotc.py [species tree file] [gene tree file]";
	print "---->\t[species tree file]: A file containing a single rooted, newick species tree.";
	print "---->\t[gene tree file]: A file containing many rooted, newick formatted gene trees, one per line. The tip labels for each gene tree must exactly match those in the species tree.\n";
	sys.exit();

specfilename = sys.argv[1];
genefilename = sys.argv[2];

if not os.path.isfile(specfilename) or not os.path.isfile(genefilename):
	sys.exit("\n** Error! Invalid file name for one of your files!\n");

sinfo, stree, sroot = tr.treeParse(open(specfilename, "r").read());
sclades = { node : set(tr.getClade(node, sinfo)) for node in sinfo if sinfo[node][2] != 'tip'};

node_counts = defaultdict(float);
lines_skipped = 0;
total_trees = 0.0;

for line in open(genefilename):
	try:
		ginfo, gtree, groot = tr.treeParse(line);
	except:
		lines_skipped += 1;
		continue;

	if len(ginfo) != len(sinfo):
		lines_skipped += 1;
		continue;

	gclades = [ set(tr.getClade(node, ginfo)) for node in ginfo if ginfo[node][2] != 'tip' ];
	for node in sclades:
		if sclades[node] in gclades:
			node_counts[node] += 1.0;
		if node == sroot and sclades[node] not in gclades:
			print len(ginfo);
			print gtree;
	total_trees += 1.0;

for node in node_counts:
	print node, round(node_counts[node]/total_trees,2);

print "-----";
print total_trees;
print lines_skipped;



