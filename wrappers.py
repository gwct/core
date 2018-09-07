#!/usr/bin/python
########################################################################################
# A script with wrapper modules for many evolutionary analysis programs.
#
# Dependencies: core
#
# Gregg Thomas, Summer 2017
########################################################################################

import sys, os, random, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core, wrapperlib as wrap

####################
parser = argparse.ArgumentParser(description="A script with wrapper modules for many evolutionary analysis programs");
parser.add_argument("-i", dest="input", help="A directory containing input files or a input file.", default=False);
parser.add_argument("-p", dest="path", help="The full path to the program you want to run. If none is specified, the script will simply try the name of the program.", default=False);
parser.add_argument("-v", dest="verbosity", help="The amount of output to print to the screen. 1 (default): print all output from the specified program, 0: print only a progress bar, program output will be redirected (maybe?).", type=int, default=1); 
parser.add_argument("-o", dest="output", help="Desired output location. If input is a file this should be a file, if a directory this will be a directory.", default=False);
# The basic options for all modules: input, path to program, verbosity, and output

parser.add_argument("--muscle", dest="muscle", help="To align the input FASTA files with MUSCLE.", action="store_true");
parser.add_argument("--pasta", dest="pasta", help="To align the input FASTA files with PASTA.", action="store_true");
parser.add_argument("--prank", dest="prank", help="To align the input FASTA files with PRANK. Used for codon alignments ONLY!", action="store_true");
parser.add_argument("--macse", dest="macse", help="To align the input FASTA files with PRANK. Used for codon alignments ONLY!", action="store_true");
parser.add_argument("--gblocks", dest="gblocks", help="To mask alignments with GBlocks.", action="store_true");
parser.add_argument("--raxml", dest="raxml", help="To make gene trees from alignments with RAxML.", action="store_true");
parser.add_argument("--codeml", dest="codeml", help="To make ancesntral reconstructions or run the branch-site test with PAML's codmel.", action="store_true");
parser.add_argument("--sdm", dest="sdm", help="To use SDM to make an average consensus distance matrix, and then R to build a neighbor-joining tree from it.", action="store_true");
parser.add_argument("--r8s", dest="r8s", help="To use r8s to estimate divergence times of a phylogeny.", action="store_true");
parser.add_argument("--notung", dest="notung", help="To use Notung to do bootstrap rearrangement or bootstrap rooting on a set of gene trees.", action="store_true");
# The various programs I have implemented wrappers for.

parser.add_argument("-seqtype", dest="seqtype", help="Specify sequence type for PASTA, GBlocks, and codeml: dna, rna, codon, or protein (default)", default="protein");
# Shared optoin between PASTA, GBlocks, and codeml.

parser.add_argument("-mode", dest="run_mode", help="GBlocks: normal (0): Accept all masked sequences. phylogenetic (1): Accept only runs that mask less than 20 percent of the original file.", type=int, default=False);
# GBlocks specific option.

parser.add_argument("-model", dest="model", help="RAxML: The DNA or protein evolution model used in RAxML runs.", default=False);
parser.add_argument("-b", dest="bootstrap_reps", help="RAxML: The number of bootstrap replicates you wish RAxML to run with its rapid bootstrapping algorithm. Default: 0", type=int, default=0);
parser.add_argument("-t", dest="threads", help="RAxML: The number of threads you wish to use for the analysis. Default: 1", type=int, default=1);
parser.add_argument("-c", dest="constraint_tree", help="RAxML: A file containing a constraint tree to be used with RAxML's -g option.", default=False);
parser.add_argument("--bl", dest="estimate_bl", help="RAxML: Use with -c to set RAxML to '-f e' to estimate branch lengths only on the constraint tree", action="store_true");
# RAxML specific options

parser.add_argument("-tree", dest="tree", help="codeml: The species tree file to use.", default=False);
parser.add_argument("--genetrees", dest="gt_opt", help="Set this to indicate that the file given by -tree contains a gene tree for each alignment to be run through codeml (in the format of the --raxml best_trees.txt file).", action="store_true");
parser.add_argument("-testbranch", dest="test_spec", help="The species you want to test with the branch-site test. Give species delimited by commas to indicate their ancestral branch.", default=False);
parser.add_argument("--prune", dest="prune", help="codeml: If not all species present in the tree will be present in each alignment, set this to prune the tree for each file. NOTE: requires Newick Utilities to be installed.", action="store_true");
parser.add_argument("--branchsite", dest="branch_site", help="codeml: By default, this program runs the null model of the branch-site test (model=2, NSsite=2, fix_omega=1, omega=1). Set this option to run the alternate model (model=2, NSsite=2, fix_omega=0, omega=1). A branch must be specified in your tree file.", action="store_true");
parser.add_argument("--anc", dest="anc", help="codeml: Set this option to perform ancestral reconstructions.", action="store_true");
parser.add_argument("--codonds", dest="codonds", help="codeml: On codon alignments calculate pairwise dn and ds rates between lineages (runmode=-2)", action="store_true");
# codeml specific options.

parser.add_argument("-numsites", dest="numsites", help="r8s: The total number of alignment columns (sites) used to build the input phylogeny.", type=int, default=False);
parser.add_argument("-calspec", dest="calspecs", help="r8s: A list of PAIRS of species that define nodes you wish to constrain times on. Species within a pair should be separated by a comma, pairs should be separated by a space (eg 'pair1s1,pair1s2 pair2s1,pair2s2').", default=False);
parser.add_argument("-calage", dest="calage", help="r8s: The calibration ages of the nodes defined by the species with -calspec. The order of this list corresponds to the order of -s. Separate ages by spaces. If constraints are to be used the keywords min and/or max are used with hyphens (eg '324 min-99.9-max-121' defines one fixed age of 324 and one constrained age).", default=False);
# r8s specific options.

parser.add_argument("-s", dest="spectree", help="Notung: A file containing a single rooted, Newick species tree for the analysis.", default=False);
parser.add_argument("--rearrange", dest="rearrange", help="Notung: Use this option to perform bootstrap rearrangement on the input gene trees.", action="store_true");
parser.add_argument("--bsroot", dest="bsroot", help="Notung: Use this option to perform bootstrap rooting on the input gene trees.", action="store_true");
parser.add_argument("-bsthresh", dest="bsthresh", help="Notung: Used with --rearrange, this specifies the bootstrap threshold for rearrangements. Default: 70", type=int, default=70);
# Notung specific options

args = parser.parse_args();
# Input option definitions.

if not args.input or not os.path.exists(args.input):
	sys.exit(core.errorOut(1, "-i must be specified and must be a valid file or directory name."));	
else:
	if os.path.isfile(args.input):
		file_flag = True;
		filelist = [os.path.abspath(args.input)];
	else:
		file_flag = False;
		filelist_init = os.listdir(args.input);
		filelist = [os.path.abspath(os.path.join(args.input, f)) for f in filelist_init];
# This checks if the input (-i) entered is valid. If so, it parses it as either a directory or a single file.

if not args.output:
	print "\n** Warning -- No output location specified. Will be determined automatically.";
if file_flag and args.verbosity == 0:
	print "\n* Message: When input type is a file, -v 0 is ignored.";
# Some messages based on the input options.

pad = 40;
if args.muscle:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --muscle, only options -i, -p, -v, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "MUSCLE", file_flag);
	print "-------------------------";
	wrap.runMuscle(filelist, file_flag, path, args.verbosity, output, logfilename);
	sys.exit();
# --muscle

if args.pasta:
	if args.seqtype not in ['dna','rna','protein']:
		sys.exit(core.errorOut(3, "PASTA accepts sequence types (-seqtype) of protein, dna, and rna"));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --pasta, only options -i, -p, -v, -seqtype, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "PASTA", file_flag);
	print core.spacedOut("Input sequence type:", pad), args.seqtype;
	print "-------------------------";
	wrap.runPasta(filelist, file_flag, path, args.seqtype, args.verbosity, output, logfilename);
	sys.exit();
# --pasta

if args.prank:
	print "\n** Warning -- PRANK is only used for codon alignments.\n"
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --prank, only options -i, -p, -v, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "PRANK", file_flag);
	print core.spacedOut("Input sequence type:", pad), args.seqtype;
	print "-------------------------";
	wrap.runPrank(filelist, file_flag, path, args.verbosity, output, logfilename);
	sys.exit();
# --prank

if args.macse:
	print "\n** Warning -- MACSE is only used for codon alignments.\n"
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --prank, only options -i, -p, -v, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "MACSE", file_flag);
	print core.spacedOut("Input sequence type:", pad), args.seqtype;
	print "-------------------------";
	wrap.runMacse(filelist, file_flag, path, args.verbosity, output, logfilename);
	sys.exit();
# --macse

if args.gblocks:
	if args.run_mode not in [0,1]:
		sys.exit(core.errorOut(4, "When running GBlocks (--gblocks) a run mode (-mode) of 0 or 1 must be specified!"));
	if args.seqtype not in ['dna','codon','protein']:
		sys.exit(core.errorOut(5, "GBlocks accepts sequence types (-seqtype) of dna, codon, and protein."));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --gblocks, only options -i, -p, -v, -seqtype, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "GBLOCKS", file_flag);
	print core.spacedOut("Input sequence type:", pad), args.seqtype;
	if args.run_mode == 0:
		print core.spacedOut("Run mode set to " + str(args.run_mode) + ":", pad), "Using default GBlocks settings (stringent).";
	elif args.run_mode == 1:
		print core.spacedOut("Run mode set to " + str(args.run_mode) + ":", pad), "Only accepting alignments with < 20% of sequence masked (for phylogenetic reconstruction).";
	print "-------------------------";
	wrap.runGblocks(filelist, file_flag, path, args.seqtype, args.run_mode, args.verbosity, output, logfilename);
	sys.exit();
# --gblocks

if args.raxml:
	if not args.model:
		sys.exit(core.errorOut(6, "When running RAxML (--raxml) a DNA or protein model (-model) must be specified!"));
	if args.bootstrap_reps < 0 or args.threads < 0:
		sys.exit(core.errorOut(7, "The number of bootstrap replicates (-b) and threads (-t) must be specified as 0 or higher!"));
	if args.constraint_tree and not os.path.isfile(args.constraint_tree):
		sys.exit(core.errorOut(8, "Invalid file name for your constraint tree (-c)."));
	if not args.constraint_tree and args.estimate_bl:
		sys.exit(core.errorOut(9, "A constraint tree (-c) must be specified with --bl"));
	# The --raxml module has several options that need to be checked for input errors before it starts.
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --raxml, only options -i, -p, -v, -model, -b, -t, -c, --bl, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "RAxML", file_flag);
	print core.spacedOut("Using DNA or protein model:", pad), args.model;
	print core.spacedOut("Number of threads:", pad), args.threads;
	if args.constraint_tree:
		args.constraint_tree = os.path.abspath(args.constraint_tree);
		print core.spacedOut("Using constraint tree:", pad), args.constraint_tree;
		if args.estimate_bl:
			print "Only estimating branch lengths on constraint tree.";
	if args.bootstrap_reps > 0:
		print "Running", args.bootstrap_reps, "bootstrap replicates on each alignment.";
	print "-------------------------";
	wrap.runRaxml(filelist, file_flag, path, args.model, args.bootstrap_reps, args.threads, args.constraint_tree, args.estimate_bl, args.verbosity, output, logfilename);
	sys.exit();
# --raxml

if args.codeml:
	if not args.path:
		sys.exit(core.errorOut(10, "A path (-p) must be specified for codeml!"));
	if args.seqtype not in ['codon','protein']:
		sys.exit(core.errorOut(11, "codeml accepts sequence types (-seqtype) of codon and protein."));
	if args.codonds and args.seqtype != 'codon':
		sys.exit(core.errorOut("11a", "The --codonds option requires -seqtype codon."))
	if not args.tree or not os.path.isfile(args.tree):
		sys.exit(core.errorOut(12, "Invalid tree file name!"));
	if args.branchsite and not args.test_spec:
		sys.ext(core.errorOut("12a", "A foreground branch must be specified to perform the branch site test. Be sure you do this for the null run as well (without --branchsite)."));
	if args.prune:
		"\n** Warning: The --prune option requires Newick Utilities to be installed and executable as 'nw_prune'!";
	# The --codeml module has several options that need to be checked for input errors before it starts.
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --codeml, only options -i, -p, -v, -seqtype, -tree, --prune, --branchsite, --anc, and -o are used!";
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "codeml", file_flag);
	print core.spacedOut("Specified sequence type:", pad), args.seqtype;
	print core.spacedOut("Using tree from file:", pad), args.tree;
	if not args.branch_site:
		print "Running NULL model for branch-site test (model=2, NSsite=2, fix_omega=1, omega=1).";
	else:
		print "Running ALTERNATE model for branch-site test model=2, NSsite=2, fix_omega=1, omega=0).";
	if args.anc:
		print "Performing ancestral reconstructions.";
	if args.prune:
		print "Pruning species tree when necessary."
	print "-------------------------";
	wrap.runCodeml(filelist, file_flag, path, args.seqtype, args.tree, args.gt_opt, args.prune, args.branch_site, args.test_spec, args.anc, args.codonds, args.verbosity, output, logfilename);
	sys.exit()
# --codeml

if args.sdm:
	if not file_flag:
		sys.exit(core.errorOut(13, "The input (-i) for --sdm must be a FILE that contains many gene trees with branch lengths (one per line, with the first line being simply the total number of trees)."));
	print "\n** Warning -- The --sdm option assumes R is callable as Rscript!\n";
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --sdm, only options -i, -p, -v, and -o are used!";
	print core.spacedOut("Using gene trees in file:", pad), args.input;
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "SDM", file_flag);
	print "-------------------------";
	wrap.runSDM(filelist, file_flag, path, args.verbosity, output, logfilename);
	sys.exit();
# --sdm

if args.r8s:
	if not file_flag:
		sys.exit(core.errorOut(14, "The input (-i) for --r8s must be a FILE that contains a single species tree with branch lengths."));
	if False in [args.numsites, args.calspecs, args.calage]:
		sys.exit(core.errorOut(15, "All of -numsites, -calspec, and -calage must be defined to run --r8s!"));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --r8s, only options -i, -p, -numsite, -calspec, -calage, and -o are used!";
	print core.spacedOut("Using species tree in file:", pad), args.input;
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "r8s", file_flag);
	wrap.runReights(filelist, file_flag, path, args.numsites, args.calspecs, args.calage, output, logfilename);
	sys.exit();
# --r8s

if args.notung:
	if file_flag:
		sys.exit(core.errorOut(-1, "Single file Notung commands not yet implemented..."));
	if (not args.rearrange and not args.bsroot) or (args.rearrange and args.bsroot):
		sys.exit(core.errorOut(16, "One and only one of --rearrange and --bsroot can be specified when using Notung."));
	if not args.spectree or not os.path.isfile(args.spectree):
		sys.exit(core.errorOut(17, "Not a valid species tree file (-s)!"));
	else:
		args.spectree = os.path.abspath(args.spectree);
	if args.bsthresh < 0 or args.bsthresh > 100:
		sys.exit(core.errorOut(18, "The bootstrap threshold (-bsthresh) for rearrangments (--rearrange) must be between 0 and 100."));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "** For --notung, only options -i, -p, -s, --rearrange, --bsroot, -v, and -o are used!";
	print core.spacedOut("Using species tree in file:", pad), args.spectree;
	path, output, logfilename = wrap.ioInfo(args.input, args.path, args.output, "notung", file_flag);
	wrap.runNotung(filelist, file_flag, path, args.spectree, args.rearrange, args.bsthresh, args.bsroot, args.verbosity, output, logfilename);
	sys.exit();
# --notung



	
