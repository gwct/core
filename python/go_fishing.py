import sys, os, scipy.stats as stats, core
import argparse

#######################################################################

def optParse(errorflag=0):
	parser = argparse.ArgumentParser(description="");
	parser.add_argument("-q", dest="query_file", help="A subset of GO terms to test for enrichment.");
	parser.add_argument("-b", dest="background_file", help="A set of GO terms used as the background for the enrichment test.");
	parser.add_argument("-a", dest="alpha_input", help="Alpha: the p-value threshold", type=float, default=0.05);
	parser.add_argument("-c", dest="correction_method", help="The multiple test correction method: Bonferroni (b), Dunn-Sidak (ds), or False discovery rate (fdr).");
	parser.add_argument("-o", dest="output_file", help="Output file to write enriched GO terms.")
	args = parser.parse_args();

	if errorflag == 0:
		if None in [args.query_file, args.background_file, args.output_file]:
			core.errorOut(1, "All of query (-q), background (-b), and output (-o) files must be specified!");
			optParse(1);

		if not os.path.exists(os.path.abspath(args.query_file)):
			core.errorOut(2, "Query file not found!");
			optParse(1);
		if not os.path.exists(os.path.abspath(args.background_file)):
			core.errorOut(3, "Background file not found!");
			optParse(1);			

		if args.alpha_input <= 0 or args.alpha_input >= 1:
			core.errorOut(4, "The p-value threshold (alpha, -a) must be between 0 and 1.");
			optParse(1);

		if args.correction_method not in ['b', 'ds', 'fdr', None]:
			core.errorOut(5, "Correction (-c) method must be one of b, ds, and fdr! Alternatively, if -c is unspecified, no correction will be done.");
			optParse(1);

		return args.query_file, args.background_file, args.alpha_input, args.correction_method, args.output_file;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

##########

def enrichedOut(go_dict, threshold):
	enriched = 0;
	for go_acc in go_dict:
		pval = go_dict[go_acc][2];
		if pval <= threshold:
			enriched += 1;
			outline = go_acc + "\t";
			for col in go_dict[go_acc]:
				outline += str(col) + "\t";
			outfile.write(outline[:-1] + "\n");

	return enriched

##########

def spacedOut(string, totlen):
#Properly adds spaces to the end of a message to make it a given length
	spaces = " " * (totlen - len(string));
	return string + spaces;

#######################################################################

infilename, gofilename, alpha, correction, outfilename = optParse();

print("# Counting total GO terms for each category...");
query_genes = [];
query_go_count = 0;
for line in open(infilename):
	if line[0] == "#":
		continue;
	line = line.strip().split("\t");
	gid = line[0];
	if gid not in query_genes:
		query_genes.append(gid);
	query_go_count += 1;


background_genes = [];
bg_go_count = 0;
for line in open(gofilename):
	if line[0] == "#":
		continue;
	line = line.strip().split("\t");
	gid = line[0];
	if gid in query_genes:
		continue;
	if gid not in background_genes:
		background_genes.append(gid);
	bg_go_count += 1;

pad = 30;

print(spacedOut("+ Total query GO terms:", pad), query_go_count, "in", len(query_genes), "genes.");
print(spacedOut("+ Total background GO terms:", pad), bg_go_count, "in", len(background_genes), "genes.");

outfile = open(outfilename, "w");
outfile.write("# GO Accession\tFamily ID\t# PS w/GO / # PS w/o GO\t# background w/ GO / # background w/o GO\t\tp-value\tOdds ratio\tGO term name\tGO domain\tGO definition\n");

print("# Running Fisher's Tests...");
done, pvals, results_dict, num_tests = [], [], {}, 0;
numlines = core.getFileLen(infilename);
i, numbars, donepercent = 0, 0, [];

for line in open(infilename):
	numbars, donepercent = core.loadingBar(i, numlines, donepercent, numbars)
	i += 1;
	if line[0] == "#":
		continue;

	line = line.strip().split("\t");
	geneid, go_acc = line[0], line[2];
	if go_acc in done:
		continue;
	done.append(go_acc);
	num_tests += 1;

	w, x, y, z = 0, 0, 0, 0;
	# w = # genes with GO term in query
	# x = # genes with GO term in background
	# y = # genes without GO term in query
	# z = # genes without GO term in background

	for next_line in open(infilename):
		if next_line[0] == "#":
			continue;
		next_line = next_line.strip().split("\t");
		next_go_acc = next_line[2];
		if next_go_acc == go_acc:
			w += 1;
	y = query_go_count - w;

	for not_line in open(gofilename):
		if not_line[0] == "#":
			continue;
		not_line = not_line.strip().split("\t");
		not_go_acc = not_line[2];
		if not_go_acc == go_acc:
			x += 1;
	z = bg_go_count - x;

	oddsratio, pvalue = stats.fisher_exact([[w,y],[x,z]], alternative='greater');
	results_dict[go_acc] = [str(w) + "/" + str(y), str(x) + "/" + str(z), pvalue, str(oddsratio), line[3], line[4], line[5]];
	pvals.append(pvalue);

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print("\n# Done!");

num_tests = float(num_tests);
print(spacedOut("+ Number of tests:", pad), int(num_tests));
print(spacedOut("+ Original alpha:", pad), alpha);
if correction == None:
	print("* Warning: Not correcting for multiple tests!");
	num_enriched = enrichedOut(results_dict, alpha);

if correction == "b":
	corrected_alpha = (alpha / num_tests);
	print(spacedOut("+ Bonferroni corrected alpha:", pad), corrected_alpha);
	num_enriched = enrichedOut(results_dict, corrected_alpha);

elif correction == "ds":
	corrected_alpha = (1 - (1 - alpha)**(1/num_tests))
	print(spacedOut("+ Dunn-Sidak corrected alpha:", pad), corrected_alpha);
	num_enriched = enrichedOut(results_dict, corrected_alpha);

elif correction == "fdr":
	pvals = sorted(pvals);
	with open("test2.txt", "w") as pfile:
		pfile.write(",".join([ str(p) for p in pvals ]));
	qvals = [];
	corrected_alpha = 0.0;
	for x in range(len(pvals)):
		pval = pvals[x];
		qval = (float(x)/num_tests) * alpha;
		qvals.append(qval);
	for x in range(len(pvals)):
		if pvals[x] < qvals[x]:
			corrected_alpha = pvals[x];

	print(spacedOut("+ FDR corrected alpha:", pad), corrected_alpha);
	num_enriched = enrichedOut(results_dict, corrected_alpha);

print(spacedOut("+ Number enriched:", pad), num_enriched);
outfile.close();
#print(pvals;