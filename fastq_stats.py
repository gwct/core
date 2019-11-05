#!/usr/bin/python3
########################################################################################
# A general purpose batch FASTQ parsing script.
#
# Dependencies: core
#
# Gregg Thomas, Fall 2019
########################################################################################

import sys, os, random, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core, fastqlib as fql

####################

def initCSV(globs):
    with open(globs['outcsv'], "w") as csvfile:
        header = ["File"];
        if globs['reads']:
            header.append("Num.reads");
        if globs['lens']:
            header.append("Avg.read.len");
        if globs['reads'] and globs['lens'] and globs['genome_size']:
            header.append("Coverage");
        if globs['bc']:
            header.extend(("A", "T", "C", "G", "N"));
        if globs['qual']:
            header.extend(("Avg.qual", "Avg.qual.first.5"));
        csvfile.write(",".join(header) + "\n");

####################
if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    parser = argparse.ArgumentParser(description="A general purpose FASTA editing script.");
    parser.add_argument("-i", dest="input", help="A directory containing FASTQ formatted files, a single FASTQ file, or a pair of FASTQ files (separated by a semicolon).");
    parser.add_argument("-g", dest="genome_size", help="The size of the original genome. If specified, coverage will be calculated.", type=int, default=False);
    parser.add_argument("-p", dest="procs", help="The number of processes the script should use. Default: 1.", type=int, default=1);
    parser.add_argument("--reads", dest="reads", help="Set to count the number of reads in each file.", default=False, action="store_true");
    parser.add_argument("--readlen", dest="lens", help="Set to calculate the average read length in each file.", default=False, action="store_true");
    parser.add_argument("--basecomp", dest="base_comp", help="Set to count the base compositions in each file.", default=False, action="store_true");
    parser.add_argument("--qual", dest="qual", help="Set to calculate the average base quality in the file and across reads.", default=False, action="store_true");
    parser.add_argument("--all", dest="all", help="Sets all of --reads, --readlen, --basecomp, --qual to True.", default=False, action="store_true");
    parser.add_argument("-outcsv", dest="outcsv", help="CSV output file for easy parsing (only log info will be printed to screen).", default=False);
    parser.add_argument("-outtxt", dest="outtxt", help="Text output file for easy reading.", default=False);
    parser.add_argument("--summary", dest="summary", help="When input is a dirctory with multiple files, set this to print only combined statistics from all files.", default=False, action="store_true");
    parser.add_argument("--overwrite", dest="overwrite", help="Set this to indicate you wish to overwrite files specified by -outcsv and -outtxt if they already exist. WARNING: This means the original contents of the file will be deleted.", default=False, action="store_true");
    #parser.add_argument("--header", dest="header", help="Set to extract the header info for each file.", default=False, action="store_true");
    args = parser.parse_args();
    # Input option definitions.

    paired = False;
    dirflag = False;

    if os.path.isdir(args.input):
        dirflag = args.input;
        infiles = [ os.path.join(args.input, f) for f in os.listdir(args.input) if any(ext in f for ext in ["fastq","fq"]) ];

        if any("_R2_" in f for f in infiles):
            paired = True;

        pairs = [];
        if paired:
            for f in infiles:
                if "R2" in f:
                    continue;
                pairs.append([f, f.replace("_R1_","_R2_")]);
            infiles = pairs;
        else:
            infiles = [[f] for f in infiles];
        # Check if input is paired.
    # If input is a directory.

    else:
        if ";" in args.input:
            paired = True;
            f = args.input.split(";");
            if not all(os.path.isfile(i) for i in f):
                sys.exit(core.errorOut(1, "Cannot find input file/directory (-i)."));

            if "R1" in f[0]:
                f1, f2 = f[0], f[1];
            else:
                f1, f2 = f[1], f[0];
            infiles = [[f1, f2]]
        # If input is paired.

        else:
            if not os.path.isfile(args.input):
                sys.exit(core.errorOut(1, "Cannot find input file/directory (-i)."));
            infiles = [[args.input]];
    # Prepare input files.

    #print(infiles);

    if args.procs < 1 or type(args.procs) != int:
        sys.exit(core.errorOut(2, "Number of processes (-p) must be an integer > 0."));
    # Make sure specified number of procs is a positive integer.

    if not args.overwrite:
        if args.outcsv and os.path.isfile(args.outcsv):
            sys.exit(core.errorOut(3, "Specified output CSV file (-outcsv) already exists. Specify a different file name or use --overwrite to overwrite the file's contents."));
        if args.outtxt and os.path.isfile(args.outtxt):
            sys.exit(core.errorOut(4, "Specified output text file (-outtxt) already exists. Specify a different file name or use --overwrite to overwrite the file's contents."));
    # Check output files.

    if args.all:
        args.reads, args.lens, args.base_comp, args.qual = True, True, True, True;

    if args.summary and not dirflag:
        print(" #" + core.getTime() + " * MESSAGE: The --summaryonly and --nosummary flags only work when input is a directory with multip fastq files. Ignoring for this run.");
        args.summary = False;

    globs = {
        'reads' : args.reads,
        'lens' : args.lens,
        'bc' : args.base_comp,
        'qual' : args.qual,
        'genome_size' : args.genome_size,
        'dirflag' : dirflag,
        'paired' : paired,
        'outcsv' : args.outcsv,
        'outtxt' : args.outtxt,
        'procs' : args.procs,
        'summary' : args.summary,
        'pyv' : sys.version[0]
    }

    print("# =======================================================================");
    print("#\t\t\t" + core.getDateTime());
    print("# Parsing FASTQ files in:\t" + args.input);

    if globs['outcsv']:
        initCSV(globs);
    if globs['outtxt']:
        with open(globs['outtxt'], "w") as txtfile:
            txtfile.write("");
    # Prepping the output files.

    fql.countReads(infiles, globs);
