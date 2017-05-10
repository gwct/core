#!/usr/bin/python
#############################################################################
#Runs RAxML on a single .fa file or a directory full of .fa files.
#
#Dependencies: core, RAxML
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
from random import randint
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Runs RAxML on a single .fa file or a directory full of .fa files. Dependencies: core, RAxML");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-r", dest="raxml_path", help="You can specify the full path to your RAxML executable here. Default: raxml (assumes you either have an alias or it is in your PATH.", default="raxml");
	parser.add_argument("-m", dest="raxml_model", help="The DNA or AA model you wish RAxML to use.");
	parser.add_argument("-b", dest="bootstrap_reps", help="The number of bootstrap replicates you wish RAxML to run with its rapid bootstrapping algorithm. Default: 0", type=int, default=0);
	parser.add_argument("-t", dest="num_threads", help="The number of threads you wish to use for the analysis. Default: 1", type=int, default=1);
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. 1: print all RAxML output, 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-c", dest="tree_combine", help="A boolean option to tell the script whether to create a file with a list of all the best trees (1) or not (0). Default: 1", type=int, default=1);
	parser.add_argument("-l", dest="log_opt", help="A boolean option to tell the script whether to create a logfile (1) or not (0). Default: 1", type=int, default=1);

	args = parser.parse_args();

	if errorflag == 0:
		if args.input == None or args.raxml_model == None:
			parser.print_help();
			sys.exit();

		if args.bootstrap_reps < 0:
			core.errorOut(1, "-b can take only positive values");
			optParse(1);

		if args.bootstrap_reps > 100:
			print " ---------------------------------------------------------------------------------------------------";
			print "|*Warning: You have specified more than 100 bootstrap replicates. This could take a very long time. |";
			print " ---------------------------------------------------------------------------------------------------";

		if args.num_threads <= 0:
			core.errorOut(2, "-t can take only positive, non-zero values");
			optParse(1);

		if args.verbosity not in [0,1]:
			core.errorOut(3, "-v must take values of either 1 or 0");
			optParse(1);

		if args.tree_combine not in [0,1]:
			core.errorOut(4, "-t must take values of either 1 or 0");
			optParse(1);

		if args.log_opt not in [0,1]:
			core.errorOut(5, "-l mus take values of either 1 or 0");
			optParse(1);

		return args.input, args.raxml_path, args.raxml_model, args.bootstrap_reps, args.num_threads, args.verbosity, args.tree_combine, args.log_opt;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

ins, rax_path, model, b, t, v, c, l = optParse(0);

starttime = core.getLogTime();

if os.path.isfile(ins):
	fileflag = 1;
	indir = os.path.dirname(os.path.realpath(ins));
	indir, script_outdir = core.getOutdir(indir, "run_raxml", starttime);
	filelist = [ins];
else:
	fileflag = 0;
	indir, script_outdir = core.getOutdir(ins, "run_raxml", starttime);
	filelist = os.listdir(indir);

bestdir = os.path.join(script_outdir, "raxml_best");
outdir = os.path.join(script_outdir, "raxml_out");

print core.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir '" + script_outdir +"'");

logfilename = os.path.join(script_outdir, "run_raxml.log");
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();

core.logCheck(l, logfilename, "=======================================================================");
core.logCheck(l, logfilename, "\t\t\tBuilding trees with RAxML");
core.logCheck(l, logfilename, "\t\t\t" + core.getDateTime());
if fileflag == 1:
	core.logCheck(l, logfilename, "INPUT    | Making tree from file:\t\t" + indir);
else:
	core.logCheck(l, logfilename, "INPUT    | Making trees from all files in:\t" + indir);
core.logCheck(l, logfilename, "INPUT    | RAxML path set to:\t\t\t" + rax_path);
core.logCheck(l, logfilename, "INFO     | Using the following DNA or AA model:\t" + model);
if b > 0:
	core.logCheck(l, logfilename, "INFO     | Performing " + str(b) + " bootstrap replicates per tree.");
else:
	core.logCheck(l, logfilename, "INFO     | Not performing bootstrap analysis.");
if t > 1:
	core.logCheck(l, logfilename, "INFO     | Using " + str(t) + " threads.");
else:
	core.logCheck(l, logfilename, "INFO     | Using 1 thread");
if v == 1:
	core.logCheck(l, logfilename, "INFO     | Printing all RAxML output to the screen.");
else:
	core.logCheck(l, logfilename, "INFO     | Silent mode. Not printing RAxML output to the screen.");
if c == 1:
	core.logCheck(l, logfilename, "INFO     | Combining RAxML best trees into a file in the input directory called raxml_best_trees.txt");
else:
	core.logCheck(l, logfilename, "INFO     | Not combining best trees to a file.");
core.logCheck(l, logfilename, "OUTPUT   | An output directory has been created within the input directory called:\t" + script_outdir);
core.logCheck(l, logfilename, "OUTPUT   | Best trees will be placed in raxml_best/, all other RAxML output will be placed in raxml_out/");
core.logCheck(l, logfilename, "-------------------------------------");

if not os.path.exists(outdir):
	core.logCheck(l, logfilename, core.getTime() + " | Creating RAxML output directory:\t" + outdir);
	cmd = "mkdir '" + outdir + "'";
	os.system(cmd);
if not os.path.exists(bestdir):
	core.logCheck(l, logfilename, core.getTime() + " | Creating RAxML best directory:\t" + bestdir);
	cmd = "mkdir '" + bestdir + "'";
	os.system(cmd);

seedFile = os.path.join(script_outdir, "raxml_seeds.txt");
core.logCheck(l, logfilename, core.getTime() + " | Creating seeds file:\t\t\t" + seedFile);
sFile = open(seedFile, "w");
sFile.write("");
sFile.close();

if b > 0:
	bseedFile = os.path.join(script_outdir, "raxml_bseeds.txt");
	core.logCheck(l, logfilename, core.getTime() + " | Creating bootstrap seeds file:\t" + bseedFile);
	sFile = open(bseedFile, "w");
	sFile.write("");
	sFile.close();

core.logCheck(l, logfilename, core.getTime() + " | Starting RAxML runs...\n");
if v == 0:
	rax_logfile = script_outdir + "raxml.log";

i = 0;
numbars = 0;
donepercent = [];

for each in filelist:
	if each.find(".fa") == -1:
		continue;

	if v == 0 and fileflag == 0:
		numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i += 1;

	if fileflag == 1:
		rax_infile = each;
		if each.find("/") != -1:
			rax_outfile = each[each.rfind("/")+1:each.index(".",each.rfind("/")+1)];
		else:
			rax_outfile = each[:each.index(".")];
	else:
		rax_infile = indir + each;
		rax_outfile = each[:each.index(".")];
	rax_outdir = os.path.join(outdir, rax_outfile + "_raxout/");

	if not os.path.exists(rax_outdir):
		os.system("mkdir '" + rax_outdir + "'");
	

	seed = str(randint(1000000,999999999));
	sFile = open(seedFile, "a");
	sline = each + "\t" + str(seed) +"\n";
	sFile.write(sline);
	sFile.close();

	if b > 0:
		boot_seed = str(randint(1000000,999999999));
		sFile = open(bseedFile, "a");
		sline = each + "\t" + str(boot_seed) +"\n";
		sFile.write(sline);
		sFile.close();
	##Generate the starting seed and bootstrap seeds (if applicable).

	rax_cmd = rax_path + " ";
	if b > 0:
		rax_cmd = rax_cmd + "-f a ";
	rax_cmd = rax_cmd + " -m " + model + " -p " + seed;
	if b > 0:
		rax_cmd = rax_cmd + " -x " + boot_seed + " -# " + str(b) + " ";
	if t > 1:
		rax_cmd = rax_cmd + "-T " + str(t) + " ";
	rax_cmd = rax_cmd + " -s '" + rax_infile + "' -n '" + rax_outfile + "' -w '" + script_outdir + "'";

	if v == 0:
		rax_cmd = rax_cmd + " >> " + rax_logfile;
	##Building the RAxML command based on the input parameters.

	if v == 1 or fileflag == 1:
		core.logCheck(l, logfilename, core.getTime() + " | RAxML Call:\t" + rax_cmd);
	else:
		lfile = open(logfilename, "a");
		lfile.write(core.getTime() + " | RAxML Call:\t" + rax_cmd + "\n");
		lfile.close();

	os.system(rax_cmd);
	##The RAxML call

	newfileList = os.listdir(script_outdir);
	for neweach in newfileList:
		if neweach.find("RAxML_bestTree") != -1:
			mv_cmd = "mv '" + script_outdir + neweach + "' '" + bestdir + "'";
			os.system(mv_cmd);
		elif neweach.find("RAxML") != -1 and neweach != "RAxML_best" and neweach != "raxml_seeds" and neweach != "RAxML_out" and neweach != "raxml_bseeds":
			mv_cmd = "mv '" + script_outdir + neweach + "' '" + rax_outdir + "'";
			os.system(mv_cmd);
	if os.path.exists(rax_infile + ".reduced"):
		mv_cmd = "mv '" + rax_infile + ".reduced '" + rax_outdir + "'";
		os.system(mv_cmd);

if v == 0:
	pstring = "100.0% complete.\n";
	sys.stderr.write('\b' * len(pstring) + pstring);


if i > 1 and c == 1:
	##Combine best trees into a single file.
	core.logCheck(l, logfilename, "\n" + core.getTime() + " | Combining best trees...");
	filelist = os.listdir(bestdir);
	tree_combine = os.path.join(script_outdir, "raxml_best_trees.txt");

	if ".DS_Store" in filelist:
		numtrees = len(filelist)-1;
	else:
		numtrees = len(filelist);

	outfile = open(tree_combine, "w");
	outfile.write(str(numtrees));
	outfile.write("\n");

	for each in filelist:
		if each == ".DS_Store":
			continue;

		infilename = os.path.join(bestdir, each);
		infile = open(infilename, "r");
		curtree = infile.read();
		infile.close();
		outfile.write(curtree);
	outfile.close();

	if b > 0:
		core.logCheck(l, logfilename, "\n" + core.getTime() + " | Combining bootstrap trees...");
		filelist = os.listdir(outdir);
		bs_filename_file = os.path.join(script_outdir, "bs_files_for_astral.txt");
		bsfile = open(bs_filename_file, "w");
		bs_tree_combine = os.path.join(script_outdir, "raxml_bs_trees.txt");
		bstreefile = open(bs_tree_combine, "w");

		for gene_dir in filelist:
			if gene_dir == ".DS_Store":
				continue;

			geneid = gene_dir[:gene_dir.index("_raxout")];
			bsfilename = os.path.join(outdir, gene_dir, "RAxML_bipartitions." + geneid);
			bstree = open(bsfilename, "r").read();
			bstreefile.write(geneid + "\t" + bstree);

			bsfile_filename = os.path.join(outdir, gene_dir, "RAxML_bootstrap." + geneid);
			bsfile.write(os.path.abspath(bsfile_filename) + "\n");

		bstreefile.close();
		bsfile.close();


core.logCheck(l, logfilename, core.getTime() + " | Done!");
core.logCheck(l, logfilename, "=======================================================================");

