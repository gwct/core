#!/usr/bin/python
#############################################################################
#Performs the likelihood ratio test for a directory full of fasta files on which
#the proper PAML runs have been made.
#
#Dependencies: core
#
#Gregg Thomas, Fall 2015
#############################################################################

import sys, os, math, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

############################################
#Function Definitions
############################################
def optParse():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", dest="input_dir", help="The directory containing your input genes.");
	parser.add_argument("-a", dest="alt_dir", help="The directory containing the PAML output from the alternate hypothesis.");
	parser.add_argument("-n", dest="null_dir", help="The directory containing the PAML output from the null hypothesis.");
	parser.add_argument("-m", dest="run_mode", help="This specifies which genes should be written to the output file: 0 = all genes, 1 = only the genes at the 1 percent significance level, 2 = only the genes at the 5 percent significance level, 3 = only the non-significant genes.", type=int, default=0);
	parser.add_argument("-o", dest="output_file", help="The prefix name of the output file. The suffix and extension (.txt) will be added based on -m.");
	args = parser.parse_args();

	if args.input_dir == None or args.alt_dir == None or args.null_dir == None or args.output_file == None:
		sys.exit(core.errorOut(1, "-i, -a, -n, and -o must all be defined"));
	if args.run_mode not in [0,1,2,3]:
		sys.exit(core.errorOut(2, "-m must take values of 1, 2, 3 or 4"));

	return args.input_dir, args.alt_dir, args.null_dir, args.run_mode, args.output_file;
#################
def getFailedFiles(logfile):
	fails = [];
	for line in open(logfile):
		if "skipped" in line:
			fails += line.strip().split(": ")[1].split(",");
	return fails;

############################################
#Main Block
############################################
#critical values: 5% = 2.71, 1% = 5.41
#conservative critical values: 5% = 3.84, 1% = 5.99
#clade test critical values: 5% = 7.82, 1% = 11.35
crit5 = 3.84;
crit1 = 5.99;
indir, altdir, nulldir, mode, outfilename = optParse();
if outfilename.endswith(".txt"):
	outfilename = outfilename.replace(".txt","");
starttime = core.getLogTime();
print "======================================================================="
print "\tPerforming likelihood ratio test (branch-site model)";
print "\t\t" + core.getDateTime();
print "INPUT    | Input directory:\t\t" + indir;
print "INPUT    | Null hypothesis directory:\t" + nulldir;
print "INPUT    | Alt hypothesis directory:\t" + altdir;
if mode == 0:
	outfilename = outfilename + "-all.txt";
	print "INFO     | Reporting results for all genes. (5% critical value = 3.84, 1% critical value = 5.99)";
elif mode == 1:
	outfilename = outfilename + "-sig1.txt";
	print "INFO     | Reporting results for only genes at the 1% critical value threshold. (5.99)";
elif mode == 2:
	outfilename = outfilename + "-sig5.txt";
	print "INFO     | Reporting results for only genes at the 5% critical value threshold. (3.84)";
elif mode == 3:
	outfilename = outfilename + "-nonsig.txt";
	print "INFO     | Reporting results for the non-significant genes.";
print "OUTPUT   | Output file:\t\t\t" + outfilename;
print "-------------------------------------"

altlogfile = os.path.join(altdir, [ f for f in os.listdir(altdir) if "run-codeml-" in f ][0] );
altfails = getFailedFiles(altlogfile);
nulllogfile = os.path.join(nulldir, [ f for f in os.listdir(nulldir) if "run-codeml-" in f ][0] );
nullfails = getFailedFiles(nulllogfile);
# Gets the files that failed to finish in the PAML runs so we can skip them.

filelist = os.listdir(indir);
# Get the list of alignment files from the input directory.

fcritcount, ocritcount, zcount, noln_genes = 0,0,0,[];
# 5 crit count, 1 crit count, zero count, genes with no likelihood score.

i, numfiles, numbars, donepercent = 0, len(filelist), 0, [];
# Variables for the progress bar.

print "Performing LRT on PAML output to test for positive selection...";
with open(outfilename, "w") as outfile:
	outfile.write("Gene ID\tnull -lnL\talt -lnL\tLR\tabove 5% critical value? (3.84)\tabove 1% critical value? (5.99)\n");
	for each in filelist:
		if each.find(".fa") == -1 or each in altfails or each in nullfails:
			continue;
		cur_file = each.replace(".fa","")
		numbars, donepercent = core.loadingBar(i, numfiles, donepercent, numbars);
		i = i + 1;

		altfilename = os.path.join(altdir, "codeml-out", cur_file + "-codemlout", cur_file + ".out");
		nullfilename = os.path.join(nulldir, "codeml-out", cur_file + "-codemlout", cur_file + ".out");

		#Reading the alt file...
		altfile = open(altfilename, "r");
		altlines = altfile.readlines();
		altfile.close();
		altlnflag = 0;
		for alt in altlines:
			if alt[:3] == "lnL":
				if alt.find("nan") != -1:
					altlnflag = 0;
				else:
					altlnL = alt[alt.index("-"):];
					altlnL = altlnL[:altlnL.index(" ")];
					altlnflag = 1;
		if altlnflag == 0:
			noln_genes.append(each);

		#Reading the null file...
		nullfile = open(nullfilename, "r");
		nulllines = nullfile.readlines();
		nullfile.close();
		nulllnflag = 0;
		for null in nulllines:
			if null[:3] == "lnL" and null.find("nan") == -1:
				if null.find("nan") != -1:
					nulllnflag = 0;
				else:
					nulllnL = null[null.index("-"):];
					nulllnL = nulllnL[:nulllnL.index(" ")];
					nulllnflag = 1;
		if nulllnflag == 0:
			if each not in noln_genes:
				noln_genes.append(each);

		if nulllnflag == 0 or altlnflag == 0:
			continue;

		negaltlnL = math.fabs(float(altlnL));
		negnulllnL = math.fabs(float(nulllnL));
		lr = 2 * (negnulllnL - negaltlnL);
		#The actual LRT

		if lr == 0:
			zcount = zcount + 1;

		fflag = False;
		oflag = False;
		outline = each + "\t" + str(negnulllnL) + "\t" + str(negaltlnL) + "\t" + str(lr);
		# Prepare the output line.

		if lr >= crit5:
			outline += "\t*";
			fcritcount += 1;
			fflag = True;
		else:
			outline += "\t";
		# Check if the LR is above the 5% critical value.

		if lr >= crit1:
			outline += "\t*";
			ocritcount += 1;
			oflag = True;
		else:
			outline += "\t";
		outline += "\n";
		# Check if the LR is above the 1% critical value.

		if mode == 0:
			outfile.write(outline);
		if mode == 1 and oflag:
			outfile.write(outline);
		if mode == 2 and fflag:
			outfile.write(outline);
		if mode == 3 and not oflag and not fflag:
			outfile.write(outline);
		# Write the line based on the results of the LRT and the current run mode.

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print "-------------------------------------"
print i, "total genes";
print fcritcount, " above 5% crit value (" + str(crit5) + ")";
print ocritcount, " above 1% crit value (" + str(crit1) + ")";
print zcount, " statistics = 0";
print float(zcount) / float(i), " percent statistics = 0";
if len(noln_genes) > 0:
	print "The following genes failed to converge in PAML:";
	for each in noln_genes:
		print each;
print "======================================================================="
