#!/usr/bin/python
#############################################################################
#Runs codeml on a single .fa file or a directory full of .fa files.
#
#Dependencies: core, treeparse, PAML, newickutils
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
from random import randint
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core
import treeparse

aas = ["G","A","V","L","I","P","F","Y","W","S","T","C","M","N","Q","K","R","H","D","E","X","-"];
nts = ["A","T","C","G","N","-","X"];

############################################
#Function Definitions
############################################
def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Runs codeml on a single .fa file or a directory full of .fa files. Dependencies: core, treeparse, PAML, newickutils");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-t", dest="tree_file", help="A user specified tree for codeml to use. If not specified, codeml will infer the tree.", default="");
	parser.add_argument("-p", dest="prune_opt", help="If not all species present in the tree will be present in each alignment, set this to 1 to prune the tree for each file. Default: 0", type=int, default=0);
	parser.add_argument("-a", dest="anc_opt", help="Option to tell PAML to do ancestral reconstruction (1) or not (0). Default: 0.", type=int, default=0);
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. 1: print all codeml output, 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-l", dest="log_opt", help="A boolean option to tell the script whether to create a logfile (1) or not (0). Default: 1", type=int, default=1);

	args = parser.parse_args();

	if args.input == None:
		parser.print_help();
		sys.exit();

	if args.prune_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 1: -p must take values of either 1 or 0 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	if args.prune_opt == 1 and args.tree_file == "":
		print " ----------------------------------------------------------";
		print "|**Error 2: With -p set to 1 a tree file must be specified |";
		print " ----------------------------------------------------------";
		parser.print_help();
		sys.exit();

	if args.anc_opt not in [0,1]:
		print " -----------------------------------------";
		print "|**Error 3: -a must take values of 1 or 0 |";
		print " -----------------------------------------";
		parser.print_help();
		sys.exit();

	if args.verbosity not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 3: -v must take values of either 1 or 0 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	if args.log_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 4: -l must take values of either 1 or 0 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	return args.input, args.tree_file, args.prune_opt, args.anc_opt, args.verbosity, args.log_opt;

#####

def logCheck(lopt, lfilename, outline):
	if lopt == 1:
		core.printWrite(lfilename, outline);
	else:
		print outline;

############################################
#Main Block
############################################

ins, treefile, prune, aopt, v, l = IO_fileParse();

starttime = core.getLogTime();
starttime = starttime.replace(":",".");

if os.path.isfile(ins):
	fileflag = 1;
	indir = os.path.dirname(os.path.realpath(ins)) + "/";
	script_outdir = ins + "-run_codeml/";
	outdir = script_outdir + "codeml_out/";
	if aopt == 1:
		ancdir = script_outdir + "anc_seqs_fa/";
	filelist = [ins];
else:
	fileflag = 0;
	filelist = os.listdir(ins);
	indir = ins;
	used = [];
	for each in filelist:
		if each.find("-run_codeml") != -1:
			used.append(int(each[:each.index("-")]));
	if used != []:
		script_outdir = ins + str(max(used)+1) + "-run_codeml-" + starttime + "/";
	else:
		script_outdir = ins + "1-run_codeml-" + starttime + "/";
	outdir = script_outdir + "codeml_out/";
	if aopt == 1:
		ancdir = script_outdir + "anc_seqs_fa/";

print core.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir " + script_outdir);

logfilename = script_outdir + "run_codeml.log";
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();

logCheck(l, logfilename, "=======================================================================");
logCheck(l, logfilename, "\t\t\tRunning codeml");
logCheck(l, logfilename, "\t\t\t" + core.getDateTime());
if fileflag == 1:
	logCheck(l, logfilename, "INPUT    | Making tree from file:\t\t" + ins);
else:
	logCheck(l, logfilename, "INPUT    | Making trees from all files in:\t" + ins);
if treefile != "":
	logCheck(l, logfilename, "INFO     | Using tree from file:\t\t" + treefile);
else:
	logCheck(l, logfilename, "INFO     | No tree file specified. codeml will infer a tree for each gene.");
if prune == 1:
	logCheck(l, logfilename, "INFO     | Pruning the tree for each gene.");
if aopt == 0:
	logCheck(l, logfilename, "INFO     | Not performing ancestral reconstruction.");
elif aopt == 1:
	logCheck(l, logfilename, "INFO     | Saving ancestral reconstructions and probabilities.");
if v == 1:
	logCheck(l, logfilename, "INFO     | Printing all codeml output to the screen.");
else:
	logCheck(l, logfilename, "INFO     | Silent mode. Not printing codeml output to the screen.");
logCheck(l, logfilename, "OUTPUT   | An output directory has been created within the input directory called:\t" + script_outdir);
logCheck(l, logfilename, "-------------------------------------");
sys.exit();
if not os.path.exists(outdir):
	logCheck(l, logfilename, core.getTime() + " | Creating codeml output directory:\t" + outdir);
	cmd = "mkdir " + outdir;
	os.system(cmd);

if aopt == 1:
	if not os.path.exists(ancdir):
		logCheck(l, logfilename, core.getTime() + " | Creating directory to pass ancestral sequences and trees:\t" + ancdir);
		cmd = "mkdir " + ancdir;
		os.system(cmd);

if prune == 1:
	logCheck(l, logfilename, core.getTime() + " | Retrieving tree info...");
	td, tree = treeparse.treeParse(open(treefile, "r").read().replace("\n",""),1);
	tips = [];
	for node in td:
		if td[node][9] == 'tip':
			tips.append(node);

logCheck(l, logfilename, core.getTime() + " | Starting codeml runs...\n");
if v == 0:
	codeml_logfile = script_outdir + "codeml.log";

ctlfilename = "codeml.ctl";

i = 0;
numbars = 0;
donepercent = [];

for each in filelist:
	if each.find(".fa") == -1:
		continue;

	if v == 0 and fileflag == 0:
		numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;

	if fileflag == 1:
		infilename = each;

	else:
		infilename = indir + each;
	gid = each[:each.index(".")];
	run_outdir = outdir + gid + "/";
	if not os.path.exists(run_outdir):
		os.system("mkdir " + run_outdir);

	if prune == 1:
		seqs = core.fastaGetDict(infilename);
		to_prune = [];
		for tip in tips:
			if (">" + tip) not in seqs:
				to_prune.append(tip);

		nw_cmd = "nw_prune " + treefile + " ";
		for tip in to_prune:
			nw_cmd = nw_cmd + tip + " ";
		nw_cmd = nw_cmd + " > pruned.tre";
#		print nw_cmd;
		os.system(nw_cmd);

	ctlFile = open(ctlfilename, "w");
	
	inline = "seqfile = " + infilename + "\n";
	ctlFile.write(inline);
	if treefile != "":
		if prune == 1:
			treeline = "treefile = pruned.tre\n";
		else:
			treeline = "treefile = " + treefile + "\n";
		ctlFile.write(treeline);
	outline = "outfile = " + run_outdir + gid + ".out\n\n";
	ctlFile.write(outline);

	ctlFile.write("noisy = 3\n");
	ctlFile.write("verbose = 0\n");
	if treefile != "":
		ctlFile.write("runmode = 0\n\n");
	else:
		ctlFile.write("runmode = 2\n\n");

	ctlFile.write("seqtype = 2\n");
	ctlFile.write("CodonFreq = 2\n");
	ctlFile.write("clock = 0\n");
	ctlFile.write("aaDist = 0\n");
	ctlFile.write("aaRatefile = /Users/Gregg/bin/paml4.7a/dat/wag.dat\n");
	ctlFile.write("model = 2\n\n");

	ctlFile.write("NSsites = 0\n\n");

	ctlFile.write("icode = 0\n");
	ctlFile.write("fix_kappa = 0\n");
	ctlFile.write("kappa = 3\n");
	ctlFile.write("fix_omega = 0\n");
	ctlFile.write("omega = 1\n\n");

	ctlFile.write("fix_alpha = 1\n");
	ctlFile.write("alpha = 0\n");
	ctlFile.write("Malpha = 0\n");
	ctlFile.write("ncatG = 10\n\n");

	ctlFile.write("getSE = 0\n");
	ctlFile.write("RateAncestor = 1\n");
	ctlFile.write("Small_Diff = .5e-6\n");

	ctlFile.close();

	codeml_cmd = "codeml " + ctlfilename;
	if v == 0:
		if os.path.isfile(ins):
			codeml_cmd = codeml_cmd + " >> " + codeml_logfile;
		else:
			codeml_cmd = codeml_cmd + " >> " + codeml_logfile;
	if v == 1 or fileflag == 1:
		logCheck(l, logfilename, core.getTime() + " | codeml Call:\t" + codeml_cmd);
	else:
		lfile = open(logfilename, "a");
		lfile.write(core.getTime() + " | codeml Call for " + infilename + ":\t" + codeml_cmd + "\n");
		lfile.close();
	os.system(codeml_cmd);

	if aopt != 0:
		anc_seqfile = ancdir + gid + "_anc.fa";
		anc_probfile = ancdir + gid + "_ancprobs.fa";
		anc_treefile = ancdir + gid + "_anc.tre";
		rstfile = open("rst", "r");
		rstlines = rstfile.readlines();
		rstfile.close();

		node_list = [];
		for k in xrange(len(rstlines)):
			if rstlines[k] == "tree with node labels for Rod Page's TreeView\n":
				anctree = rstlines[k+1].replace(" ","");
				atfile = open(anc_treefile, "w");
				atfile.write(anctree);
				atfile.close();

			if rstlines[k] == "List of extant and reconstructed sequences\n":
				asfile = open(anc_seqfile, "w");

				j = 4;
				while rstlines[k+j].find("Overall accuracy of the") == -1:
					if rstlines[k+j] != "\n":
						tmpline = rstlines[k+j].replace("\n", "");
						tmpline = tmpline.split("  ");
						tmpline = filter(None, tmpline);
						node_list.append(tmpline[0]);

						title = ">" + tmpline[0] + "\n";
						asfile.write(title);
						seq = tmpline[1].replace(" ","") + "\n";
						asfile.write(seq);
					j = j + 1;
				if aopt == 1:
					asfile.close();
				break;		

		curseqs = {};
		for n in node_list:
			curseqs[n] = "";

		for k in xrange(len(rstlines)):
			if rstlines[k] == "Prob of best state at each node, listed by site\n":
				j = 4;
				while rstlines[k+j] != "\n":
					tmpline = filter(None, rstlines[k+j].replace("\n","").split("  "));
					cur_states = tmpline[2].split(": ");
					extant = list(cur_states[0].replace(" ",""));
					ancs = filter(None, cur_states[1].split(" "));
					all_states = extant + ancs;

					for n in xrange(len(node_list)):
						cur_spec = node_list[n];
						cur_aa = all_states[n];
						if cur_aa.find("(") != -1 and cur_aa.find(")") != -1:
							cur_aa = cur_aa.replace("(","_").replace(")","") + " ";
						curseqs[cur_spec] = curseqs[cur_spec] + cur_aa;
					j = j + 1;

		asfile = open(anc_probfile, "w");
		for seq in curseqs:
			title = ">" + seq + "\n";
			asfile.write(title);
			for base in curseqs[seq]:
				asfile.write(base);
			asfile.write("\n");
		asfile.close();

	
	newfileList = os.listdir(os.getcwd());
	for neweach in newfileList:
		if neweach in ["2NG.dN","2NG.dS","2NG.t","codeml.ctl","lnf","rst","rst1","rub","pruned.tre"]:
			mv_cmd = "mv " + neweach + " " + run_outdir;
			os.system(mv_cmd);


if v == 0 and fileflag == 0:
	pstring = "100.0% complete.\n";
	sys.stderr.write('\b' * len(pstring) + pstring);
logCheck(l, logfilename, core.getTime() + " | Done!");
logCheck(l, logfilename, "=======================================================================");




