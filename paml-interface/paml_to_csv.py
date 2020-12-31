#!/usr/bin/python
############################################################
# Reads results from codeml output M2 output per branch
############################################################

import sys, os, argparse, lib.pamlcore as pcore

############################################################
# Options

parser = argparse.ArgumentParser(description="Parse PAML codeml output for a tip and target branches");
parser.add_argument("-i", dest="input", help="Directory containing subdirectories of codeml runs.", default=False);
parser.add_argument("-m", dest="model", help="The PAML model that was used to generate the files in -i. Default: m1", default="m1");
parser.add_argument("-d", dest="meta", help="A file containing metadata. Tab delimited with columns: id, feature type, chromosome, start coord, end coord, strand. Directories in -i must be formatted <id>-*", default=False);
parser.add_argument("-o", dest="output", help="An output .csv file.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit(" * Error 1: Please provide a valid input directory (-i).");

if args.model not in ["m1", "m2", "cmc"]:
    sys.exit(" * Error 2: Model (-m) must be one of: m1, m2, cmc");

if args.meta and not os.path.isfile(args.meta):
    sys.exit(" * Error 3: Cannot find meta data file: " + args.meta);

if not args.output:
    sys.exit(" * Error 2: Please provide the name of an output file (-o).")

if os.path.isfile(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output file (-o) already exists! Explicity specify --overwrite to overwrite it.");

pad = 25;

with open(args.output, "w") as outfile:
    pcore.runTime("# codeml m1 output parser", outfile);
    pcore.PWS("# IO OPTIONS", outfile);
    pcore.PWS(pcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    pcore.PWS(pcore.spacedOut("# PAML model:", pad) + args.model, outfile);
    if args.meta:
        pcore.PWS(pcore.spacedOut("# Metadata file:", pad) + args.meta, outfile);
    pcore.PWS(pcore.spacedOut("# Output file:", pad) + args.output, outfile);
    if args.overwrite:
        pcore.PWS(pcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous output file.", outfile);
    pcore.PWS("# ----------------", outfile);

    features = False;
    if args.meta:
        pcore.PWS("# " + pcore.getDateTime() + " Reading metadata file: " + args.meta, outfile);
        features = pcore.readMeta(args.meta);
        pcore.PWS(pcore.spacedOut("# Features read:        ", pad) + str(len(features)), outfile);               
        pcore.PWS("# ----------------", outfile);
    # Read the feature metadata.

    pcore.PWS("# " + pcore.getDateTime() + " Begin parsing PAML output...", outfile);
    if args.model == "m1":
        import lib.m1 as m1;
        m1.parse(args.input, features, outfile, pad);

    if args.model == "cmc":
        import lib.cmc as cmc;
        cmc.parse(args.input, features, outfile, pad);

    if args.model == "m2":
        import lib.m2 as m2;
        m2.parse(args.input, features, outfile, pad);
    # Load the library for the model used and pass everything to it.


    

