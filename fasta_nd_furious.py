#!/usr/bin/python
########################################################################################
# A general purpose batch FASTA editing script.
#
# Dependencies: core
#
# Gregg Thomas, Summer 2017
########################################################################################

import sys, os, random, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core, fastalib as fa

####################
parser = argparse.ArgumentParser(description="A general purpose FASTA editing script.");
parser.add_argument("-i", dest="input", help="A directory containing FASTA formatted files or a single FASTA file.");
parser.add_argument("--countpos", dest="count_pos", help="With this flag, the total number of positions in all files will be counted.", default=False, action="store_true");
parser.add_argument("--countaln", dest="count_aln", help="With this flag, the total number of columns in an alignment will be counted, along with other alignment stats.", default=False, action="store_true");
parser.add_argument("--concat", dest="concat", help="Set this option to concatenate all alignments in the input directory.", default=False, action="store_true");
parser.add_argument("--combine", dest="combine", help="Set this option to combine all fasta sequences in a directory. A single file may contain one or more sequences.", default=False, action="store_true");
parser.add_argument("--split", dest="split", help="Given an input FASTA file, this option will split all sequences into individual files. Be sure to define -delim otherwise the entire header will be used as the output file name!", default=False, action="store_true");
parser.add_argument("--trim", dest="trim", help="Use this to trim the FASTA headers at the first occurrence of a character defined by -delim", default=False, action="store_true");
parser.add_argument("-extract", dest="extract", help="Given a file of FASTA titles, this will extract those titles and sequences from the FASTA file specified with -i.", default=False);
parser.add_argument("-relabel", dest="relabel", help="Use this option to place a new label in the FASTA headers. New labels are specified with -newlabel. 1 (default): place at beginning of old header, 2: completely replace old header, 3: place at end of old header. If -delim is specified, that will be used as the seperator. The default delimiter is _", type=int, default=False);
parser.add_argument("-delim", dest="header_delim", help="A character at which the FASTA headers will be trimmed. For a delimter of ' ' (space), please enter 'space'", default=False);
parser.add_argument("-newlabel", dest="new_label", help="The new label(s) to be added to the FASTA headers when -relabel is specified. If a single label is specified, this will be used for all sequences. If a dict like list of labels is specified (ie 'oldlabel1:newlabel1,oldlabel2:newlabel2') then the old headers will be searched for the old labels and new labels will be placed appropriately.", default=False);
parser.add_argument("-rmseq", dest="label_rm", help="A comma delimited list of labels. Any sequence with any of those labels in their header will be removed.", default=False);
parser.add_argument("-rmstart", dest="start_rm", help="This will remove the start M or ATG from input sequences. Enter sequence type: 'protein' or 'codon'", default=False);
parser.add_argument("-replace", dest="replace", help="This option will replace all characters in each sequence with another character. For example, AB will replace all As with Bs. If the input is an alignment and A: is entered, all As will be replaced with another base/aa that is not present in the column. For multiple replacements, enter as: AB,CD,EF", default=False);
parser.add_argument("-outfile", dest="outfile", help="The output file for commands: --concat, --combine, -extract");
parser.add_argument("-outdir", dest="outdir", help="The output directory for commands:");
args = parser.parse_args();
# Input option definitions.

if args.header_delim == 'space':
	args.header_delim = " ";
# The special case of the -delim option for a space character.

if args.input == None or not os.path.exists(args.input):
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

if all(opt == False for opt in [args.count_pos, args.count_aln, args.concat, args.combine, args.split, args.trim, args.relabel, args.label_rm, args.start_rm, args.replace, args.extract]):
	print "\n** Warning: No options specified! Just running --countpos.";
	args.count_pos = True;

if args.count_pos:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "Counting the total number of positions (AAs or NTs) in:\t", args.input;
	fa.countPos(filelist);
	sys.exit();
# --countpos

if args.count_aln:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "Counting stats from all alignments in:\t", args.input;
	fa.countAln(filelist);
	sys.exit();
# --countaln

if args.concat or args.combine:
	if file_flag:
		sys.exit(core.errorOut(2, "-i must be a directory with multiple FASTA files when used with --concat or --combine."));
if args.concat or args.combine or args.extract:
	if args.outfile == None:
		sys.exit(core.errorOut(3, "-outfile must be set when using the --concat or --combine or -extract options."));
# This check makes sure the proper input and output have been specified for the --concat and --combine options:
# input must be a directory with multiple .fa files, output must be a single file (-outfile).

if args.concat:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Concatenating alignments in:", 40), args.input;
	print core.spacedOut("Writing concatenated alignments to:", 40), args.outfile;
	if args.header_delim == False:
		print "\n** Warning: No FASTA header delimiter set (-delim). Full FASTA headers will be used for concatenation.";
	fa.concat(filelist, args.header_delim, args.outfile);
	sys.exit();
# --concat

if args.combine:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Combining all .fa files in:", 40), args.input;
	print core.spacedOut("Writing combined sequences to:", 40), args.outfile + "\n";
	fa.combine(filelist, args.outfile);
	sys.exit();
# --combine

if args.split:
	if not file_flag:
		sys.exit(core.errorOut(4, "-i must be a single FASTA file with multiple sequences in it when used with --split"));
	if not args.outdir:
		sys.exit(core.errorOut(5, "-outdir must be set when using the --split option"));
	if not args.header_delim:
		print "\n** Warning: You have not specified a FASTA header delimiter to determine output file names. File names will be as long as your headers!\n";
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Splitting all sequences in:", 40), args.input;
	print core.spacedOut("Delimiting headers for file names at:", 40), args.header_delim;
	print core.spacedOut("Writing individual files to:", 40), args.outdir;
	fa.split(filelist, args.header_delim, args.outdir);
	sys.exit();
## --split

if any([args.trim, args.relabel, args.label_rm, args.start_rm, args.replace]) and file_flag and args.outdir:
	sys.exit(core.errorOut(4, "You have specified an input file with an output directory (-outdir). Please specify and output file (-outfile) instead!"));
elif any([args.trim, args.relabel, args.label_rm, args.start_rm, args.replace]) and not file_flag and args.outfile:
	sys.exit(core.errorOut(5, "You have specified an input directory with an output file (-outfile). Please specify and output directory (-outdir) instead!"));
else:
	if file_flag:
		output = args.outfile;
	else:
		output = args.outdir;
		if output != None and not os.path.isdir(output):
			print "\n++ Creating output directory...";
			os.system("mkdir " + output);
	if output == None:
		output = args.input;
		print "\n** Warning: You have not specified an output destination! Your input location/name will be used (with an added label) for all output files!\n"
# For the rest of the options, users can input either a file or a directory. The output option must match (-outfile or -outdir). This ensures
# that it does. If no output option is entered, the output destination is determined based on the input file name.

if args.trim:
	if args.header_delim == None:
		sys.exit(core.errorOut(7, "-delim must be set with --trim"));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Trimming all .fa files in:", 45), args.input;
	print core.spacedOut("Writing trimmed (.trim.fa) sequences to:", 45), output;
	fa.trim(filelist, args.header_delim, file_flag, output);
	sys.exit();
## --trim

if args.relabel:
	if args.relabel not in [1,2,3]:
		sys.exit(core.errorOut(8, "-relabel must take values of 1, 2, or 3"));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Relabeling all .fa files in:", 50), args.input;
	print core.spacedOut("Writing relabeled (.relabel.fa) sequences to:", 50), output;
	fa.relabel(filelist, args.relabel, args.header_delim, args.new_label, file_flag, output);
	sys.exit();
## -relabel

if args.label_rm:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Removing sequences from all .fa files in:", 60), args.input;
	print core.spacedOut("Removing all sequences with the following labels:", 60), args.label_rm;
	print core.spacedOut("Writing filtered (.rmlab.fa) sequences to:", 60), output;
	fa.removeSeq(filelist, args.label_rm, file_flag, output);
	sys.exit();
## -rmseq

if args.start_rm:
	if args.start_rm.lower() not in ['protein','prot','p','codon','c']:
		sys.exit(core.errorOut(10, "Sequence types for removing starts (-rmstart) must be 'protein' or 'codon'"));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Removing all starts from:", 50), args.input;
	print core.spacedOut("Writing all sequences (.rmstart) to:", 50), output;
	fa.removeStarts(filelist, args.start_rm.lower(), file_flag, output);
	sys.exit();	
## -rmstart

if args.replace:
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Doing replacements of bases in:", 50), args.input;
	print core.spacedOut("Doing the following replacements:", 50), args.replace;
	print core.spacedOut("Writing replaced (.repl) sequences to:", 50), output;
	fa.replaceBase(filelist, args.replace, file_flag, output);
	sys.exit();
# -replace

if args.extract:
	print args.extract;
	print os.path.isfile(args.extract);
	if os.path.isfile(args.extract):
		titles = open(args.extract, "r").read().split("\n");
	else:
		titles = args.extract.split(",");
	print titles;
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Extracting sequences in:", 50), args.input;
	print core.spacedOut("Extracting the following titles: ", 50) + ",".join(titles);
	print core.spacedOut("Writing extracted sequences to:", 50), output;
	fa.extractSeqs(filelist, titles, args.header_delim, file_flag, output);
	sys.exit();
# -extract





