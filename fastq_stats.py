#!/usr/bin/python
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
parser = argparse.ArgumentParser(description="A general purpose FASTA editing script.");
parser.add_argument("-i", dest="input", help="A directory containing FASTQ formatted files, a single FASTQ file, or a pair of FASTQ files (separated by a semicolon).");
parser.add_argument("-g", dest="genome_size", help="The size of the original genome. If specified, coverage will be calculated.", type=int, default=False);
parser.add_argument("--reads", dest="reads", help="Set to count the number of reads in each file.", default=False, action="store_true");
parser.add_argument("--readlen", dest="lens", help="Set to calculate the average read length in each file.", default=False, action="store_true");
parser.add_argument("--basecomp", dest="base_comp", help="Set to count the base compositions in each file.", default=False, action="store_true");
parser.add_argument("--qual", dest="qual", help="Set to calculate the average base quality in the file and across reads.", default=False, action="store_true");
parser.add_argument("--header", dest="header", help="Set to extract the header info for each file.", default=False, action="store_true");
parser.add_argument("-outfile", dest="outfile", help="The output file for commands: --concat, --combine, -extract");
parser.add_argument("-outdir", dest="outdir", help="The output directory for commands:");
args = parser.parse_args();
# Input option definitions.

paired = False;
dirflag = False;

if os.path.isdir(args.input):
    dirflag = True;
    infiles = [ os.path.join(args.input, f) for f in os.listdir(args.input) if any(ext in f for ext in ["fastq","fq"]) ];

    if any("R2" in f for f in infiles):
        paired = True;

    pairs = [];
    if paired:
        for f in infiles:
            if "R2" in f:
                continue;
            pairs.append([f, f.replace("R1","R2")]);
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

#print(infiles);

globs = {
    'reads' : args.reads,
    'lens' : args.lens,
    'bc' : args.base_comp,
    'qual' : args.qual,
    'genome_size' : args.genome_size,
    'dirflag' : dirflag,
    'paired' : paired,
    'pyv' : sys.version[0]
}

print("=======================================================================");
print("\t\t\t" + core.getDateTime());
print("Parsing FASTQ files in:\t" + args.input);
fql.countReads(infiles, globs);
