#!/usr/bin/python
############################################################
# For generators, when logcheck returns some runs with errors
# this can be used to comment out the files that completed 
# successfully.
# Inputs: 1. A file with the log file names of scripts that
#            errors (logcheck | cut -d ":" -f 1 | uniq > file)
#         2. The original job script generated.
############################################################

import sys, os, argparse, core, re, gzip
from datetime import datetime

############################################################

def unComment(f):
    outlines = [];
    begin_cmds = False;
    for line in open(f):
        if begin_cmds:
            if line[0] == "#":
                line = line[1:];
            if line[0] == " ":
                line = line[1:];

        if line == "# BEGIN CMDS\n":
            begin_cmds = True;

        outlines.append(line);

    with open(f, "w") as outfile:
        for outline in outlines:
            outfile.write(outline);

############################################################
parser = argparse.ArgumentParser(description="Job script modifications");
parser.add_argument("-j", dest="jobfile", help="The job script you wish to modify.", default=False);
parser.add_argument("-m", dest="modfile", help="The file containing log files for the lines you want to modify. All other lines will be commented out.", default=False);
parser.add_argument("-r", dest="repl", help="A modification you wish to make to the commands. Format: 'old cmd,new cmd'", default=False);
parser.add_argument("--anti", dest="anti", help="Set to comment all lines IN the modfile (-m)", action="store_true", default=False);
parser.add_argument("--uncomment", dest="uncomment", help="Set this option to simply uncomment every line in the job file.", action="store_true", default=False);
args = parser.parse_args();

core.runTime("# Comments out and/or modifies lines in a job bash script");

if not args.jobfile:
    sys.exit(" * ERROR: Please specify a job file to modify with -j.");
if not os.path.isfile(args.jobfile):
    sys.exit(" * ERROR: Job file not found: " + args.jobfile);
# Check for the job file.

if args.uncomment:
    print("--uncomment set: simply uncommenting all the commands in jobfile: " + args.jobfile);
    unComment(args.jobfile);
    print("# Done!");
    print("# ----------------");
    sys.exit();
# The --uncomment option.

if not args.modfile:
    sys.exit(" * ERROR: Please specify a mod file to modify with -m.");
if not os.path.isfile(args.modfile):
    sys.exit(" * ERROR: Mod file not found: " + args.modfile);
# Check for the mod file.

print("# Job file: " + args.jobfile);
print("# Mod file: " + args.modfile);

if args.repl:
    if ", " in args.repl:
        args.repl = args.repl.replace(", ", ",");
    repl = args.repl.split(",");
    print("# Replacing: " + repl[0]);
    print("# With:      " + repl[1]);
# Parse the replace option.

modlines = list(filter(None, open(args.modfile, "r").read().split("\n")));
#print(modlines);

outlines = [];
begin_cmds = False;
for line in open(args.jobfile):
    if begin_cmds:
        if not args.anti:
            comment_current = True;
            for modline in modlines:
                if modline in line:
                    comment_current = False;

                    if args.repl and repl[0] in line:
                        line = line.replace(repl[0], repl[1]);

        elif args.anti:
            comment_current = False;
            for modline in modlines:
                if modline in line:
                    comment_current = True;

                    if args.repl and repl[0] in line:
                        line = line.replace(repl[0], repl[1]);

        if comment_current and line[0] != "#":
            line = "# " + line;

    if line == "# BEGIN CMDS\n":
        begin_cmds = True;

    outlines.append(line);

with open(args.jobfile, "w") as outfile:
    for outline in outlines:
        outfile.write(outline);

print("# Done!");
print("# ----------------");