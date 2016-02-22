#!/usr/bin/python
#############################################################################
#Performs several steps of the supertree building process utilizing many programs
#
#Dependencies: core, SDM, R with the ape package, newickutils (for rooting trees),
#r8s (for smoothing trees)
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
from subprocess import Popen, PIPE
import core

#"Apis_mellifera,Bombus_impatiens Pachypsylla_venusta,Frankliniella_occidentallis Tribolium_castaneum,Athalia_rosae Anopheles_gambiae,Lutzomyia_longipalpis Catajapyx_aquilonaris,Centruroides_sculpturatus Athalia_rosae,Pachypsylla_venusta Hapegnathos_saltator,Dufourea_novaeangliae Ladona_fulva,Ephemera_danica Frankliniella_occidentallis,Pediculus_humanus"
#25.7,309.15,228,99.9,117.55,228,91.85,228,222

#"Ladona_fulva,Blatella_germanica Frankliniella_occidentalis,Pachypsylla_venusta Athalia_rosae,Atta_cephalotes Athalia_rosae,Lucilia_cuprina Plutella_xylostella,Heliconius_melpomene Anopheles_gambiae,Drosophila_melanogaster Pachypsylla_venusta,Halyomorpha_halys Apis_mellifera,Atta_cephalotes Apis_mellifera,Bombus_terrestris"
#324,222,228,309.15,124.2,99.9,99.9,91.85,25.7

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_file", help="A file containing a list of trees on which to run SDM and NJ OR a file containing a single tree on which to run r8s.");
	parser.add_argument("-r", dest="r_output_file", help="A file name for R to write the Neighbor Joining tree.");

	parser.add_argument("-j", dest="nj_opt", help="A boolean option to use SDM to create a consensus matrix and R to create a NJ tree. Default: 0.", type=int, default=0);
	parser.add_argument("-o", dest="nj_outgroup", help="The outgroup by which the NJ tree will be rooted.");
	parser.add_argument("-t", dest="reroot_opt", help="Boolean to reroot (1) the NJ tree or not (0). If set to 1, -o must also be specified. Default: 0", type=int, default=0);

	parser.add_argument("-d", dest="div_est_opt", help="A boolean option to estimate divergence times from the NJ tree with r8s (1) or not (0). Default: 0.", type=int, default=0);
	parser.add_argument("-e", dest="r8s_output_file", help="A file name for r8s to write the final output.");
	parser.add_argument("-n", dest="num_sites", help="The total number of sites from the alignments used to make the tree; used by r8s.");
	parser.add_argument("-s", dest="cal_specs", help="A list of PAIRS of species that define nodes you wish to constrain times on. Species within a pair should be separated by a comma, pairs should be separated by a space (eg 'pair1s1,pair1s2 pair2s1,pair2s2').");
	parser.add_argument("-a", dest="cal_age", help="The calibration ages of the nodes defined by the species in -s. The order of this list corresponds to the order of -s. Separate ages by commas. If constraints are to be used the keywords min and/or max are used with hyphens (eg '324,min-99.9-max-121' defines one fixed age of 324 and one constrained age).");
	parser.add_argument("-l", dest="log_opt", help="A boolean option to tell the script whether to create a logfile (1) or not (0). Default: 1", type=int, default=1);

	parser.add_argument("-z", dest="script_logdir", help="A directory in which to place the script output directory. If none is specified, this will default to the directory of the input file");
	parser.add_argument("-x", dest="logdir_suffix", help="A string to add on to the end of the output directory.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.input_file == None or args.nj_opt == None or args.nj_opt not in [0,1] or args.div_est_opt == None or args.div_est_opt not in [0,1]:
			core.errorOut(1, "-i must always be defined. One of -j or -d must also always be defined as 1");
			optParse(1);

		if args.reroot_opt not in [0,1]:
			core.errorOut(2, "-t must take values of either 1 or 0");
			optParse(1);

		if args.reroot_opt == 1 and args.nj_outgroup == None:
			core.errorOut(3, "-When -t is set to 1, an outgroup must be specified with -o");
			optParse(1);

		if args.div_est_opt not in [0,1]:
			core.errorOut(4, "-d must take values of either 1 or 0");
			optParse(1);

		elif args.div_est_opt == 1:
			if args.r8s_output_file == None or args.num_sites == None or args.cal_specs == None or args.cal_age == None:
				core.errorOut(5, "You are missing one or more of the options for div time estimation with r8s. -e, -n, -s, and -a must all be defined");
				optParse(1);
			else:
				if args.cal_specs.find(" ") != -1 and args.cal_age.find(",") != -1:
					cal_specs = args.cal_specs.split(" ");
					cal_age = args.cal_age.split(",");
				else:
					cal_specs = [args.cal_specs];
					cal_age = [args.cal_age];
				if len(cal_specs) != len(cal_age):
					core.errorOut(6, "You must enter the same number of calibration nodes (-s) and calibration ages (-a)");
					optParse(1);

		else:
			args.r8s_output_file = None;
			args.num_sites = None;
			cal_specs = None;
			args.cal_age = None;

		if args.log_opt not in [0,1]:
			core.errorOut(7, "-l must take values of either 1 or 0");
			optParse(1);

		return args.input_file, args.r_output_file, args.nj_opt, args.nj_outgroup, args.reroot_opt, args.div_est_opt, args.r8s_output_file, args.num_sites, cal_specs, cal_age, args.log_opt, args.script_logdir, args.logdir_suffix;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

infilename, routfilename, njopt, outgroup, rr, d, r8soutfilename, numsites, calspec, calage, l, script_outdir_initial, outdir_suffix = optParse(0);

starttime = core.getLogTime();

if infilename.find("/") != -1:
	indir = os.path.dirname(os.path.realpath(infilename)) + "/";
	infilename = infilename[infilename.rfind("/")+1:];
else:
	indir = os.getcwd() + "/";

indir, script_outdir = core.getOutdir(indir, "supertreemaker", starttime);
print script_outdir;
print os.path.basename(os.path.normpath(script_outdir));
if script_outdir_initial != None:
	if not os.path.isdir(script_outdir_initial):
		core.errorOut(8, "-z must be a valid directory");
		optParse(1);

	script_outdir = os.path.join(script_outdir_initial, os.path.basename(os.path.normpath(script_outdir)));
if outdir_suffix != None:
	if script_outdir[-1] == "/":
		script_outdir = script_outdir[:len(script_outdir)-1] + "-" + outdir_suffix + "/";
	else:
		script_outdir = script_outdir + "-" + outdir_suffix + "/";

print core.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir '" + script_outdir + "'");

logfilename = script_outdir + "supertreemaker.log";
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();
##Pre-run prep: creating log files and output directories, etc...

core.logCheck(l, logfilename, "=======================================================================");
core.logCheck(l, logfilename, "\tSupertree making with SDM, R, newickutils, and r8s");
core.logCheck(l, logfilename, "\t\t\t" + core.getDateTime());
core.logCheck(l, logfilename, "INPUT    | Making tree from file:\t\t\t" + infilename);
core.logCheck(l, logfilename, "INPUT    | Input file located in:\t\t\t" + indir);
if njopt == 1:
	core.logCheck(l, logfilename, "INFO     | Using Average Consensus method in SDM to build distance matrix.");
	core.logCheck(l, logfilename, "INFO     | Using R to build a neighbor-joining tree from matrix.");
else:
	core.logCheck(l, logfilename, "INFO     | Not creating consensus tree.");
if rr == 1:
	core.logCheck(l, logfilename, "INFO     | Rooting the NJ tree with species:\t\t" + outgroup);
if d == 0:
	core.logCheck(l, logfilename, "INFO     | NOT estimating divergence times with r8s.");
else:
	core.logCheck(l, logfilename, "INFO     | Will estimate divergence times with r8s.");
	core.logCheck(l, logfilename, "INFO     | Number of sites from alignments:\t\t" + numsites);
	core.logCheck(l, logfilename, "INFO     | Calibrating with the following nodes and ages:\n---------");
	core.logCheck(l, logfilename, core.spacedOut("NODE 1", 40) + core.spacedOut("NODE 2", 40) + "AGE");
	for n in xrange(len(calspec)):
		calspec[n] = calspec[n].split(",");
		core.logCheck(l, logfilename, core.spacedOut(calspec[n][0], 40) + core.spacedOut(calspec[n][1], 40) + calage[n]);
	core.logCheck(l, logfilename, "---------");
core.logCheck(l, logfilename, "OUTPUT   | An output directory has been created:\t" + script_outdir);
if njopt == 1:
	core.logCheck(l, logfilename, "OUTPUT   | Writing NJ tree to:\t\t\t\t" +  routfilename);
if d == 1:
	core.logCheck(l, logfilename, "OUTPUT   | Writing r8s output to:\t\t\t" + r8soutfilename);
core.logCheck(l, logfilename, "-------------------------------------");
##Info block

if njopt == 1:
	sdmfilename = script_outdir + infilename + "_sdm_mat.txt";
	tmpfilename = script_outdir + "/rmat.tmp";
	routfilename = script_outdir + routfilename;
	##Some more file prep

	core.logCheck(l, logfilename, core.getTime() + " | Running Average Consensus within SDM...");
	sdm_cmd = "java -jar ~/bin/SDM/SDM.jar -i " + indir + infilename + " -d ACS97 -t T -f Phylip_square";
	print sdm_cmd;
	os.system(sdm_cmd);
	##Runs SDM with ACS97

	core.logCheck(l, logfilename, core.getTime() + " | Moving SDM output files to main output directory...");
	mv_cmd = "mv " + indir + "*sdm* " + script_outdir;
	os.system(mv_cmd);
	##Moves all SDM output files to script outdir

	core.logCheck(l, logfilename, core.getTime() + " | Reading and parsing distance matrix...");
	infile = open(sdmfilename, "r");
	inlines = infile.readlines();

	specs = [];

	for x in xrange(len(inlines)):
		if x == 0 or inlines[x] == "\n":
			continue;
		while inlines[x].find("  ") != -1:
			inlines[x] = inlines[x].replace("  "," ");
		specs.append(inlines[x][:inlines[x].index(" ")]);

	specs = " ".join(specs);

	tmpfile = open(tmpfilename,"w");
	tmpfile.write(specs);
	tmpfile.write("\n");
	for x in xrange(len(inlines)):
		if x == 0 or inlines[x] == "\n":
			continue;
		tmpfile.write(inlines[x]);
	tmpfile.close();
	##Re-writes the SDM distance matrix to be readable by R

	core.logCheck(l, logfilename, core.getTime() + " | Calling R to make NJ tree...\n");
	rcmd = "Rscript ~/bin/core/corelib/nj_tree.r " + tmpfilename + " " + routfilename;
	os.system(rcmd);
	##Runs the NJ algorithm within the ape package of R

	routfile = open(routfilename, "r");
	nj_tree = routfile.read().replace("\n","");
	routfile.close();
	core.logCheck(l, logfilename, "\n ----Unrooted NJ tree----");
	core.logCheck(l, logfilename, nj_tree);
	##Reads and prints the resulting tree
else:
	nj_tree = open(indir + infilename, "r").read().replace("\n","");

if rr == 1:
##Uses newickutils to root the tree at the outgroup if specified
	core.logCheck(l, logfilename, "\n" + core.getTime() + " | Rooting NJ tree with newickutils...\n");
	rootedtreefile = routfilename[:routfilename.index(".tre")] + "_rooted.tre";
	nwcmd = "nw_reroot " + routfilename + " " + outgroup + " > " + rootedtreefile;
	os.system(nwcmd);
	nj_tree = open(rootedtreefile,"r").read().replace("\n","");
	core.logCheck(l, logfilename, "\n ----Rooted NJ tree----");
	core.logCheck(l, logfilename, nj_tree);


if d == 1:
##Uses r8s to estimate divergence times on the rooted tree
	core.logCheck(l, logfilename, "\n" + core.getTime() + " | Preparing tree and input file for r8s...");
	#r8soutfilename = script_outdir + r8soutfilename;
	r8soutfilename = os.path.join(script_outdir, r8soutfilename);

	#r8sinputfile = script_outdir + "r8s_input_file.txt";
	r8sinputfile = os.path.join(script_outdir, "r8s_input_file.txt");
#	nwcmd = "nw_prune nwtmp.txt " + outgroup + " > nwtmp.txt";
#	p = Popen(nwcmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
#	nj_tree, stderr = p.communicate()
#	nj_tree = nj_tree.replace("\n","");

	#calnode = calspec[0][:3] + calspec[1][:3];
	calnodes = [];
	for node in calspec:
		calnodes.append(node[0][:3] + node[1][:3]);

	r8sfile = open(r8sinputfile, "w");
	r8sfile.write("#NEXUS\nbegin trees;\n");
	outline = "tree nj_tree = [&U] " + nj_tree + "\nEnd;\n";
	r8sfile.write(outline);
	r8sfile.write("begin rates;\n");
	outline = "blformat nsites=" + numsites + " lengths=persite ultrametric=no;\n";
	r8sfile.write(outline);
	r8sfile.write("collapse;\n");
	for n in xrange(len(calnodes)):
		outline = "mrca " + calnodes[n] + " " + calspec[n][0] + " " + calspec[n][1] + ";\n";
		r8sfile.write(outline);
	for n in xrange(len(calnodes)):
		if calage[n].find("min") != -1 and calage[n].find("max") != -1:
			curcal = calage[n].split("-");
			outline = "constrain taxon=" + calnodes[n] + " minage=" + curcal[1] + " maxage=" + curcal[3] +";\n";
		elif calage[n].find("min") != -1:
			curcal = calage[n].split("-");
			outline = "constrain taxon=" + calnodes[n] + " minage=" + curcal[1] +";\n";
		elif calage[n].find("max") != -1:
			curcal = calage[n].split("-");
			outline = "constrain taxon=" + calnodes[n] + " minage=" + curcal[1] +";\n";
		else:
			outline = "fixage taxon=" + calnodes[n] + " age=" + calage[n] + ";\n";
		r8sfile.write(outline);
	r8sfile.write("divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;\n");
	r8sfile.write("describe plot=chronogram;\n");
	r8sfile.write("describe plot=tree_description;\n");
	r8sfile.write("end;");

	r8sfile.close();

	core.logCheck(l, logfilename, core.getTime() + " | Calling r8s to smooth the tree...");
	r8scmd = "r8s -b -f '" + r8sinputfile + "' > '" + r8soutfilename + "'";
	os.system(r8scmd);

	for line in open(r8soutfilename):
		continue;
	div_tree = line;
	for node in calnodes:
		div_tree = div_tree.replace(node,"");
	div_tree = div_tree[div_tree.index("= ")+2:];
	core.logCheck(l, logfilename, "\n ----Smoothed tree----");
	core.logCheck(l, logfilename, div_tree);

core.logCheck(l, logfilename, core.getTime() + " | Done!");
core.logCheck(l, logfilename, "=======================================================================");
