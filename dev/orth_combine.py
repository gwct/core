#!/usr/bin/python
########################################################################################
# Part of my ENSEMBL CAFE pipeline. This takes two (or more?) files containing 
# lists of orthologous genes and their confidences (either 1 or 0), selects the
# one-to-one orthologs, and combines them. The first column in each file should
# be the same species. Necessary because Ensembl only allows you to export orthologs
# from 6 species at a time, so when working with more than 7 species, two or more
# lists are required.
#
# Input file format:
# Main species ID	Species 1 ID	Species 1 confidence	Species 2 ID	Species 2 confidence	 .... etc.
#
# Gregg Thomas, Spring 2015
########################################################################################

import sys, os, argparse, itertools, core

####################

print "\n###### Hello ######"
print "orth_combine call: " + " ".join(sys.argv) + "\n"

parser = argparse.ArgumentParser();
parser.add_argument("-i", dest="input", help="A comma delimited LIST of Ensembl ortholog lists to combine. The first column of each file must be the same species.", default=False);
parser.add_argument("-o", dest="output", help="Output file name for combined ortholog list.", default=False);
args = parser.parse_args();

if args.input:
	infiles = args.input.split(",");
	for infile in infiles:
		if not os.path.isfile(infile):
			sys.exit(" ** ERROR 1: input file not found: ", infile);
else:
	sys.exit(" ** ERROR 2: at least one input file must be given (-1)");
# Check the input files

if not args.output:
	sys.exit(" ** ERROR 3: output file not speciefied (-o)");
## I/O parsing

print "# =======================================================================";
print "# Combining files:\t\t\t", infiles;
print "# Writing output to:\t\t\t" + args.output;
print "# Writing only possible one-to-one orthologs to output file.";
print "# -------------------------------------";
print "# " + core.getTime() + " Reading files...";

file_lines = {};
key_ids = {};

for infile in infiles:
	file_lines[infile] = [];
	key_ids[infile] = [];

	i = 0;
	for line in open(infile):
		if i == 0:
			line = line.split("\t");
			numcols = len(line);
			numspec = (numcols - 1) / 2;
			i = i + 1;
			continue;

		line = line.strip().split("\t");
		key_id = line[0];
		if key_id == '' or len(line) != numcols:
			continue;

		cur_conf = [];
		cur_specs = [];
		for j in range(len(line)):
			if j != 0 and j % 2 == 0:
				cur_conf.append(line[j]);
			else:
				cur_specs.append(line[j]);

		if cur_specs.count('') == 0 and cur_conf.count('1') == len(cur_conf):
			file_lines[infile].append(cur_specs);
			key_ids[infile].append(key_id);



print "# -------------------------------------";
print "# " + core.getTime() + " Retrieving key IDs in both files...";

comb_ids = list(set.intersection(*map(set,key_ids.values())));
num_ids = len(comb_ids);

print "# -------------------------------------";
print "# " + core.getTime() + " Retrieving lines with shared key IDs (all combinations)...";

i = 0;
numbars = 0;
donepercent = [];
written = [];
with open(args.output, "w") as outfile:
	for cid in comb_ids:
		numbars, donepercent = core.loadingBar(i, num_ids, donepercent, numbars);
		i += 1;
		cur_id = {};
		for infile in infiles:
			cur_id[infile] = [];
			for line in file_lines[infile]:
				if line[0] == cid:
					cur_id[infile].append(line);

		skip = False;
		for infile in infiles:
			if len(cur_id[infile]) > 1 and cur_id[infile].count(cur_id[infile][0]) != len(cur_id[infile]):
					skip = True; 
		if skip:
			#print cur_id;
			continue;
		# This skips if the query ID has multiple lines in a file that do not match, meaning one or more species has
		# two orthologs with this query.

		id_vals = [];
		for each in cur_id:
			id_vals.append(cur_id[each]);

		for b in itertools.product(*id_vals):
			outline = b[0] + b[1][1:];
			outline = "\t".join(outline) + "\n";
			if outline not in written:
				outfile.write(outline);
				written.append(outline);

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\n# " + core.getTime() + " Done!";
print "# =======================================================================";

