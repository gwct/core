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
import sys, argparse, itertools
import core

############################################
#Function Definitions
############################################

def optParse(errorflag):
# This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_list", help="A comma delimited LIST of Ensembl ortholog lists to combine. The first column of each file must be the same species.");
	parser.add_argument("-o", dest="output_file", help="Output file name for combined ortholog list.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.input_list == None or args.output_file == None:
			core.errorOut(1, "Both -i and -o must be defined");
			optParse(1);

		return args.input_list.split(","), args.output_file;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

infiles, outfilename = optParse(0);

#print infiles, outfilename;

print "# =======================================================================";
print "# Combining files:\t\t\t", infiles;
print "# Writing output to:\t\t\t" + outfilename;
print "# Writing only possible one-to-one orthologs to output file.";
print "# -------------------------------------";
print "# " + core.getTime() + " Reading files...";

file_lines = {};
key_ids = {};

for each in infiles:
	infile = open(each, "r");
	inlines = infile.readlines();
	infile.close();

	file_lines[each] = [];
	key_ids[each] = [];

	i = 0;
	for line in inlines:
		if i == 0:
			line = line.split("\t");
			numcols = len(line);
			numspec = (numcols - 1) / 2;

			i = i + 1;
			continue;

		line = line.replace("\n","").split("\t");
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
			file_lines[each].append(cur_specs);
			key_ids[each].append(key_id);

print "# -------------------------------------";
print "# " + core.getTime() + " Retrieving key IDs in both files...";

comb_ids = list(set.intersection(*map(set,key_ids.values())));
num_ids = len(comb_ids);

print "# -------------------------------------";
print "# " + core.getTime() + " Retrieving lines with shared key IDs (all combinations)...";

i = 0;
numbars = 0;
donepercent = [];

outfile = open(outfilename, "w");
for cid in comb_ids:
	numbars, donepercent = core.loadingBar(i, num_ids, donepercent, numbars);
	i += 1;
	cur_id = {};
	for each in infiles:
		cur_id[each] = [];
		for line in file_lines[each]:
			if line[0] == cid:
				cur_id[each].append(line);

	id_vals = [];
	for each in cur_id:
		id_vals.append(cur_id[each]);

	for b in itertools.product(*id_vals):
		outline = b[0] + b[1][1:];
		outline = "\t".join(outline) + "\n";
		outfile.write(outline);

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\n# " + core.getTime() + " Done!";
print "# =======================================================================";

