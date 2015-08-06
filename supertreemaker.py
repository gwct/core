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

############################################
#Function Definitions
############################################

def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_file", help="A file containing a list of trees onw which to run SDM and NJ.");
	parser.add_argument("-r", dest="r_output_file", help="A file name for R to write the Neighbor Joining tree.");

	parser.add_argument("-o", dest="nj_outgroup", help="The outgroup by which the NJ tree will be rooted.");
	parser.add_argument("-t", dest="reroot_opt", help="Boolean to reroot (1) the NJ tree or not (0). If set to 1, -o must also be specified. Default: 0", type=int, default=0);

	parser.add_argument("-d", dest="div_est_opt", help="A boolean option to estimate divergence times from the NJ tree with r8s (1) or not (0). Default: 0.", type=int, default=0);
	parser.add_argument("-e", dest="r8s_output_file", help="A file name for r8s to write the final output.");
	parser.add_argument("-n", dest="num_sites", help="The total number of sites from the alignments used to make the tree; used by r8s.");
	parser.add_argument("-s", dest="cal_specs", help="A comma delimited (no spaces) list of the two species in the NJ tree you wish to calibrate ages by.");
	parser.add_argument("-a", dest="cal_age", help="The calibration age of the node defined by the species in -s.");
	parser.add_argument("-l", dest="log_opt", help="A boolean option to tell the script whether to create a logfile (1) or not (0). Default: 1", type=int, default=1);

	args = parser.parse_args();

	if args.input_file == None or args.r_output_file == None:
		parser.print_help();
		sys.exit();

	if args.reroot_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 1: -t must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		sys.exit();

	if args.reroot_opt == 1 and args.nj_outgroup == None:
		print " -----------------------------------------------------------------------";
		print "|**Error 2: When -t is set to 1, an outgroup must be specified with -o. |";
		print " -----------------------------------------------------------------------";

	if args.div_est_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 3: -d must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	elif args.div_est_opt == 1:
		if args.nj_outgroup == None or args.r8s_output_file == None or args.num_sites == None or args.cal_specs == None or args.cal_age == None:
			print " ----------------------------------------------------------------------------------------------------------------------";
			print "|**Error 4: You are missing some options for div time estimation with r8s. -o, -e, -n, -s, and -a all must be defined. |";
			print " ----------------------------------------------------------------------------------------------------------------------";
			parser.print_help();
			sys.exit();
		else:
			cal_specs = args.cal_specs.split(",");
			args.reroot_opt = 1;

	else:
		args.r8s_output_file = None;
		args.num_sites = None;
		cal_specs = None;
		args.cal_age = None;

	if args.log_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 5: -l must take values of either 1 or 0 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();


	return args.input_file, args.r_output_file, args.nj_outgroup, args.reroot_opt, args.div_est_opt, args.r8s_output_file, args.num_sites, cal_specs, args.cal_age, args.log_opt;

#####

def logCheck(lopt, lfilename, outline):
	if lopt == 1:
		core.printWrite(lfilename, outline);
	else:
		print outline;

############################################
#Main Block
############################################

infilename, routfilename, outgroup, rr, d, r8soutfilename, numsites, calspec, calage, l = IO_fileParse();

starttime = core.getLogTime();

if infilename.find("/") != -1:
	indir = os.path.dirname(os.path.realpath(infilename)) + "/";
	infilename = infilename[infilename.rfind("/")+1:];
else:
	indir = os.getcwd() + "/";

script_outdir = indir + starttime + "-supertreemaker/";

print core.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir " + script_outdir);

logfilename = script_outdir + "supertreemaker.log";
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();
##Pre-run prep: creating log files and output directories, etc...

logCheck(l, logfilename, "=======================================================================");
logCheck(l, logfilename, "\tSupertree making with SDM, R, newickutils, and r8s");
logCheck(l, logfilename, "\t\t\t" + core.getDateTime());
logCheck(l, logfilename, "INPUT    | Making tree from file:\t\t\t" + infilename);
logCheck(l, logfilename, "INPUT    | Input file located in:\t\t\t" + indir);
logCheck(l, logfilename, "INFO     | Using Average Consensus method in SDM to build distance matrix.");
logCheck(l, logfilename, "INFO     | Using R to build a neighbor-joining tree from matrix.");
if rr == 1:
	logCheck(l, logfilename, "INFO     | Rooting the NJ tree with species:\t\t" + outgroup);
if d == 0:
	logCheck(l, logfilename, "INFO     | NOT estimating divergence times with r8s.");
else:
	logCheck(l, logfilename, "INFO     | Will estimate divergence times with r8s.");
	logCheck(l, logfilename, "INFO     | Number of sites from alignments:\t\t" + numsites);
	logCheck(l, logfilename, "INFO     | Calibrating at the node defined by:\t\t" + str(calspec));
	logCheck(l, logfilename, "INFO     | Calibration age at that node:\t\t" + calage);
	logCheck(l, logfilename, "INFO     | Writing r8s output to:\t\t\t" + r8soutfilename);
logCheck(l, logfilename, "OUTPUT   | An output directory has been created:\t" + script_outdir);
logCheck(l, logfilename, "OUTPUT   | Writing NJ tree to:\t\t\t\t" +  routfilename);
if d == 1:
	logCheck(l, logfilename, "OUTPUT   | Writing r8s output to:\t\t\t" + r8soutfilename);
logCheck(l, logfilename, "-------------------------------------");
##Info block

sdmfilename = script_outdir + infilename + "_sdm_mat.txt";
tmpfilename = script_outdir + "/rmat.tmp";
routfilename = script_outdir + routfilename;
if r8soutfilename != None:
	r8soutfilename = script_outdir + r8soutfilename;
##Some more file prep

logCheck(l, logfilename, core.getTime() + " | Running Average Consensus within SDM...");
sdm_cmd = "java -jar ~/bin/SDM/SDM.jar -i " + indir + infilename + " -d ACS97 -t T -f Phylip_square";
print sdm_cmd;
os.system(sdm_cmd);
##Runs SDM with ACS97

logCheck(l, logfilename, core.getTime() + " | Moving SDM output files to main output directory...");
mv_cmd = "mv " + indir + "*sdm* " + script_outdir;
os.system(mv_cmd);
##Moves all SDM output files to script outdir

logCheck(l, logfilename, core.getTime() + " | Reading and parsing distance matrix...");
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

logCheck(l, logfilename, core.getTime() + " | Calling R to make NJ tree...\n");
rcmd = "Rscript ~/bin/core/corelib/nj_tree.r " + tmpfilename + " " + routfilename;
os.system(rcmd);
##Runs the NJ algorithm within the ape package of R

routfile = open(routfilename, "r");
nj_tree = routfile.read().replace("\n","");
routfile.close();
logCheck(l, logfilename, "\n ----Unrooted NJ tree----");
logCheck(l, logfilename, nj_tree);
##Reads and prints the resulting tree

if rr == 1:
##Uses newickutils to root the tree at the outgroup if specified
	logCheck(l, logfilename, "\n" + core.getTime() + " | Rooting NJ tree with newickutils...\n");
	rootedtreefile = routfilename[:routfilename.index(".tre")] + "_rooted.tre";
	nwcmd = "nw_reroot " + routfilename + " " + outgroup + " > " + rootedtreefile;
	os.system(nwcmd);
	nj_tree = open(rootedtreefile,"r").read().replace("\n","");
	logCheck(l, logfilename, "\n ----Rooted NJ tree----");
	logCheck(l, logfilename, nj_tree);
	

if d == 1:
##Uses r8s to estimate divergence times on the rooted tree
	logCheck(l, logfilename, "\n" + core.getTime() + " | Preparing tree and input file for r8s...");

	r8sinputfile = script_outdir + "r8s_input_file.txt";
#	nwcmd = "nw_prune nwtmp.txt " + outgroup + " > nwtmp.txt";
#	p = Popen(nwcmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
#	nj_tree, stderr = p.communicate()
#	nj_tree = nj_tree.replace("\n","");

	calnode = calspec[0][:3] + calspec[1][:3];

	r8sfile = open(r8sinputfile, "w");
	r8sfile.write("#NEXUS\nbegin trees;\n");
	outline = "tree nj_tree = [&U] " + nj_tree + "\nEnd;\n";
	r8sfile.write(outline);
	r8sfile.write("begin rates;\n");
	outline = "blformat nsites=" + numsites + " lengths=persite ultrametric=no;\n";
	r8sfile.write(outline);
	r8sfile.write("collapse;\n");
	outline = "mrca " + calnode + " " + calspec[0] + " " + calspec[1] + ";\n";
	r8sfile.write(outline);
	outline = "fixage taxon=" + calnode + " age=" + calage + ";\n";
	r8sfile.write(outline);
	r8sfile.write("divtime method=pl algorithm=qnewt cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;\n");
	r8sfile.write("describe plot=chronogram;\n");
	r8sfile.write("describe plot=tree_description;\n");
	r8sfile.write("end;");

	r8sfile.close();

	logCheck(l, logfilename, core.getTime() + " | Calling r8s to smooth the tree...");
	r8scmd = "r8s -b -f " + r8sinputfile + " > " + r8soutfilename;
	os.system(r8scmd);

	for line in open(r8soutfilename):
		continue;
	div_tree = line[line.index("("):].replace(calnode,"");
	logCheck(l, logfilename, "\n ----Smoothed tree----");
	logCheck(l, logfilename, div_tree);

logCheck(l, logfilename, core.getTime() + " | Done!");
print "=======================================================================";
