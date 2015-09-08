#!/usr/bin/python
#############################################################################
#Runs GBlocks on a single .fa file or a directory full of .fa files.
#
#Dependencies: core, GBlocks
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Runs GBlocks on a single .fa file or a directory full of .fa files. Dependencies: core, GBlocks");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-r", dest="gblocks_path", help="You can specify the full path to your GBlocks executable here. Default: gblocks (assumes you either have an alias or it is in your PATH.", default="gblocks");
	parser.add_argument("-t", dest="seq_type", help="Choose from: protein (p, default), dna, (d), or codon (c).", default="p");
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. 1: print all GBlocks output, 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-l", dest="log_opt", help="A boolean option to tell the script whether to create a logfile (1) or not (0). Default: 1", type=int, default=1);

	args = parser.parse_args();

	if errorflag == 0:
		if args.input == None:
			parser.print_help();
			sys.exit();

		st = args.seq_type.lower();
		if st not in ["p","d","c","protein","dna","codon"]:
			core.errorOut(1, "-t must take values of p, d, or c");
			optParse(1);

		if len(st) > 1:
			st = st[:1];

		if args.verbosity not in [0,1]:
			core.errorOut(2, "-v must take values of either 1 or 0");
			optParse(1);

		if args.log_opt not in [0,1]:
			core.errorOut(3, "-l must take values of either 1 or 0");
			optParse(1);

		return args.input, args.gblocks_path, st, args.verbosity, args.log_opt;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

#####

def logCheck(lopt, lfilename, outline):
	if lopt == 1:
		core.printWrite(lfilename, outline);
	else:
		print outline;

############################################
#Main Block
############################################

ins, gb_path, seqtype, v, l = optParse(0);

starttime = core.getLogTime();

if os.path.isfile(ins):
	fileflag = 1;
	indir = os.path.dirname(os.path.realpath(ins)) + "/";
	outdir = indir + starttime + "-gblocks/";
	filelist = [ins];
else:
	fileflag = 0;
	indir, outdir = core.getOutdir(ins, "gblocks", starttime);
	filelist = os.listdir(indir);

print core.getTime() + " | Creating main output directory...";
os.system("mkdir " + outdir);

logfilename = outdir + "run_gblocks.log";
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();

logCheck(l, logfilename, "=======================================================================");
logCheck(l, logfilename, "\t\t\tMasking alignments with GBlocks");
logCheck(l, logfilename, "\t\t\t" + core.getDateTime());
if fileflag == 1:
	logCheck(l, logfilename, "INPUT    | Masking alignment from file: " + ins);
else:
	logCheck(l, logfilename, "INPUT    | Masking alignments from all files in: " + indir);
logCheck(l, logfilename, "INFO     | GBlocks path set to: " + gb_path);
logCheck(l, logfilename, "INFO     | Sequence type set to: " + seqtype);
if v == 1:
	logCheck(l, logfilename, "INFO     | Printing all GBlocks output to the screen.");
else:
	logCheck(l, logfilename, "INFO     | Silent mode. Not printing GBlocks output to the screen.");
logCheck(l, logfilename, "OUTPUT   | An output directory has been created within the input directory called: " + outdir);
logCheck(l, logfilename, "-------------------------------------");

num_aligns = 0;
for a in filelist:
	if a.find(".fa") != -1:
		num_aligns = num_aligns + 1;

if v == 0:
	gb_logfile = outdir + "gblocks.log";

logCheck(l, logfilename,  core.getTime() + " | Runnning GBlocks on " + str(num_aligns) + " alignments...");

acc_mask = 0;

i = 0;
numbars = 0;
donepercent = [];
for each in filelist:
	if v == 0 and fileflag == 0:
		numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;

	if each.find(".fa") == -1:
		continue;

	if fileflag == 1:
		infilename = each;
		if each.find("/") != -1:
			gb_outfile = each[each.rfind("/")+1:each.index(".fa")] + "-gb.fa";
		else:
			gb_outfile = each[:each.index(".fa")] + "-gb.fa";
	else:
		infilename = indir + each;
		gb_outfile = each[:each.index(".fa")] + "-gb.fa";

	inseqs = core.fastaGetDict(infilename);
	seqlen = len(inseqs[inseqs.keys()[0]]);
	b1 = int(round(0.5 * len(inseqs))) + 1;

	gb_cmd = "gblocks " + infilename + " -t=" + seqtype + " -b1=" + str(b1) + " -b2=" + str(b1) + " -b3=" + str(seqlen) + " -b4=2 -b5=a";
	if v == 0:
		gb_cmd = gb_cmd + " >> " + gb_logfile;

	if v == 1 or fileflag == 1:
		logCheck(l, logfilename, core.getTime() + " | GBlocks Call:\t" + gb_cmd);
	else:
		lfile = open(logfilename, "a");
		lfile.write(core.getTime() + " | GBlocks Call:\t" + gb_cmd + "\n");
		lfile.close();
	os.system(gb_cmd);

	gbseqs = core.fastaGetDict(infilename+"-gb");
	gblen = len(gbseqs[gbseqs.keys()[0]]);
	difflen = seqlen - gblen;
	percdiff = float(difflen) / float(seqlen) * 100.0;
	percgb = float(gblen) / float(seqlen) * 100.0;


	outline = core.getTime() + " | Percent of original sequence masked:\t" + str(percdiff);

	if percdiff <= 20:
		outline = outline + ". Accepting masked alignment.";
		acc_mask = acc_mask + 1;
		mv_cmd = "mv " + infilename + "-gb " + outdir + gb_outfile;
		os.system(mv_cmd);

	else:
		outline = outline + ". Rejecting masked alignment.";
		cp_cmd = "cp " + infilename + " " + outdir + gb_outfile;
		os.system(cp_cmd);

		rm_cmd = "rm " + infilename + "-gb";
		os.system(rm_cmd);

	if v == 1 or fileflag == 1:
		logCheck(l, logfilename, outline);
	else:
		lfile = open(logfilename, "a");
		lfile.write(outline + "\n");
		lfile.close();

	rm_cmd = "rm " + infilename + "-gb.htm ";
	os.system(rm_cmd);

if v == 0 and fileflag == 0:
	pstring = "100.0% complete.\n";
	sys.stderr.write('\b' * len(pstring) + pstring);
logCheck(l, logfilename, core.getTime() + " | Done!");
logCheck(l, logfilename, core.getTime() + " | " + str(acc_mask) + " out of " + str(num_aligns) + " masked alignments accepted.");
logCheck(l, logfilename, "=======================================================================");
