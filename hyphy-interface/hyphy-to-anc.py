#!/usr/bin/python
############################################################
# Reads results from Hyphy json output
############################################################

import sys, os, argparse, lib.hpcore as hpcore

############################################################
# Options

parser = argparse.ArgumentParser(description="Parse Hyphy json output");
parser.add_argument("-i", dest="input", help="Directory containing hyphy json output files with ancestral reconstructions (currently from SLAC).", default=False);
parser.add_argument("-m", dest="model", help="The Hyphy model that was used to generate the files in -i. Default: slac", default="slac");
#parser.add_argument("-d", dest="meta", help="A file containing metadata. Tab delimited with columns: id, feature type, chromosome, start coord, end coord, strand. Directories in -i must be formatted <id>-*", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory. Will be created if it doesn't exist.", default=False);
parser.add_argument("-g", dest="outgroups", help="A comma separated list of OUTGROUP SPECIES (tips). Will be used to root gene trees.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit(" * Error 1: Please provide a valid input directory (-i).");

if args.model not in ["slac"]:
    sys.exit(" * Error 2: Model (-m) must be one of: slac");

# if args.meta and not os.path.isfile(args.meta):
#     sys.exit(" * Error 3: Cannot find meta data file: " + args.meta);

if not args.outgroups:
    sys.exit(" * Error 2: Please provide at least one outgroup species to root gene trees (-g).");
while ", " in args.outgroups:
    args.outgroups = args.outgroups.replace(", ", ",");
outgroups = args.outgroups.split(",");

if not args.output:
    sys.exit(" * Error 3: Please provide the name of an output directory (-o).");

if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 4: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
if not os.path.isdir(args.output):
    os.system("mkdir " + args.output);

base_out = os.path.basename(os.path.normpath(args.output));
#output_file = os.path.join(globs['outdir'], base_out + ".txt");
logfilename = os.path.join(args.output, base_out + ".log");
# Main output and log file

pad = 25;

with open(logfilename, "w") as logfile:
    hpcore.runTime("# HyPhy output parser", logfile);
    hpcore.PWS("# IO OPTIONS", logfile);
    hpcore.PWS(hpcore.spacedOut("# Input directory:", pad) + args.input, logfile);
    hpcore.PWS(hpcore.spacedOut("# Hyphy model:", pad) + args.model, logfile);
    # if args.meta:
    #     hpcore.PWS(hpcore.spacedOut("# Metadata file:", pad) + args.meta, logfile);
    hpcore.PWS(hpcore.spacedOut("# Outgroups:", pad) + args.outgroups, logfile);
    hpcore.PWS(hpcore.spacedOut("# Output file:", pad) + args.output, logfile);
    if args.overwrite:
        hpcore.PWS(hpcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous output file.", logfile);
    hpcore.PWS("# ----------------", logfile);

    # features = False;
    # if args.meta:
    #     hpcore.PWS("# " + hpcore.getDateTime() + " Reading metadata file: " + args.meta, logfile);
    #     features = hpcore.readMeta(args.meta);
    #     hpcore.PWS(hpcore.spacedOut("# Features read:        ", pad) + str(len(features)), logfile);               
    #     hpcore.PWS("# ----------------", logfile);
    # Read the feature metadata.

    hpcore.PWS("# " + hpcore.getDateTime() + " Begin parsing Hyphy output...", logfile);
    if args.model == "slac":
        import lib.slac as slac;
        slac.parseAnc(args.input, args.output, outgroups, logfile, pad);
    # Load the library for the model used and pass everything to it.
