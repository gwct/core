#############################################################################
# Parsing and printing the options and meta-info.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import lib.core as CORE

#############################################################################

def inputPathCheck(path, path_type, req, opt, globs, glob_key):
# A function to check input files and directories

    if path_type == "file":
        checker = os.path.isfile;
    elif path_type == "dir":
        checker = os.path.isdir;
    # Get the right type of path to check

    if req and not path:
        CORE.errorOut("OPC1", opt + " must be provided.", globs);
    # If the path is required and not provided, error out here

    if path:
        if not checker(path):
            CORE.errorOut("OPC1", "Cannot find " + opt + ": " + path, globs);
        # If the path doesn't exist, error out here

        else:
            globs[glob_key] = os.path.abspath(path);
        # Otherwise add the provided path to the global params dict

    return globs;

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py

    try:
        import psutil
        globs['psutil'] = True;
    except:
        globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    parser = argparse.ArgumentParser(description="Generates concatenated dS, dN, and dN/dS estimates across all genes for a branch");
    parser.add_argument("-i", dest="tree_info_file", help="csv file with species tree info, including clades for each node.", default=False);
    parser.add_argument("-r", dest="rates_dir", help="The directory of csv files from a hyphy SLAC run.", default=False);
    parser.add_argument("-f", dest="filter_file", help="The file with the loci to excldue from the analysis.", default=False);
    parser.add_argument("-s", dest="subset_file", help="Subset of genes to include in the analysis.", default=False);
    parser.add_argument("-o", dest="output_dir", help="The name of the directory to output files.", default=False);
    parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
    # I/O params and options

    parser.add_argument("-n", dest="num_procs", help="The number of processes to use. Default: 1", default=False);
    # User params

    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    #parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent isofilter from reporting detailed information about each step.", action="store_true", default=False);
    # Run options

    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    args = parser.parse_args();

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    globs['overwrite'] = args.overwrite;
    # Check run mode options.

    globs = inputPathCheck(args.tree_info_file, "file", True, "-i", globs, "tree-info-file");
    # Check the tree info file

    globs = inputPathCheck(args.rates_dir, "dir", True, "-r", globs, "csv-rate-dir");
    # Check the rates dir

    globs = inputPathCheck(args.filter_file, "file", False, "-f", globs, "filter-file");
    globs = inputPathCheck(args.subset_file, "file", False, "-s", globs, "subset-file");
    # Check the optional filer and subset files

    if not args.output_dir:
        CORE.errorOut("OP6", "Must provide an output directory name defined with -o.", globs);
    elif args.output_dir and os.path.isdir(args.output_dir) and not args.overwrite:
        CORE.errorOut("OP7", "Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.", globs);
    elif args.output_dir and not os.path.isdir(args.output_dir):
        os.system("mkdir " + args.output_dir);
    globs['outdir'] = args.output_dir;
    # Output directory

    if args.num_procs:
        if not CORE.isPosInt(args.num_procs):
            CORE.errorOut("OP10", "The number of processes defined with -n must be a positive integer.", globs);
        globs['num-procs'] = int(args.num_procs);
    # The number of processes to use

    globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
    globs['output-file'] = os.path.join(globs['outdir'], globs['run-name'] + ".csv");
    globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
    # Main output and log file

    if args.quiet_flag:
        globs['quiet'] = True;
    # Check the quiet option

    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    startProg(globs);
    # After all the essential options have been set, call the welcome function.

    return globs;

#############################################################################

def startProg(globs):
# A nice way to start the program.
    print("#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Welcome to branch_avgs -- Averaging rates across a species tree while accounting for discordance.");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Version " + globs['version'] + " released on " + globs['releasedate']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# branch_avgs was developed by " + globs['authors']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + CORE.getDateTime());
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 30;
    opt_pad = 30;
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Tree info file:", pad) + globs['tree-info-file']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Rates directory:", pad) + globs['csv-rate-dir']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Output directory:", pad) + globs['outdir']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Output file:", pad) + globs['output-file']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));
    # Input/Output
    #######################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Option", pad) + CORE.spacedOut("Current setting", opt_pad) + "Current action");      

    if globs['filter-file']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -f", pad) +
                    CORE.spacedOut(globs['filter-file'], opt_pad) + 
                    "Loci listed in this file will be excluded from analysis");
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -f", pad) +
                    CORE.spacedOut("NA", opt_pad) + 
                    "No filter file provided - no loci will be filtered from analysis");         
    # Reporting the optional filter file   

    if globs['subset-file']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -s", pad) +
                    CORE.spacedOut(globs['subset-file'], opt_pad) + 
                    "Only loci listed in this file will be included in analysis");
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -s", pad) +
                    CORE.spacedOut("NA", opt_pad) + 
                    "No subset file provided - no loci will be subset for analysis");    
    # Reporting the optional subset file

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -p", pad) +
                CORE.spacedOut(str(globs['num-procs']), opt_pad) + 
                "branch_avgs will use this many processes in parallel");  
    # Reporting the number of processes

    if globs['overwrite']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --overwrite", pad) +
                    CORE.spacedOut("True", opt_pad) + 
                    "branch_avgs will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option

    if not globs['quiet']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("False", opt_pad) + 
                    "Time, memory, and status info will be printed to the screen while isofilter is running.");
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "No further information will be printed to the screen while isofilter is running.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
    # Reporting the quiet option

    # if globs['debug']:
    #     CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --debug", pad) + 
    #                 CORE.spacedOut("True", opt_pad) + 
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option.

    if globs['norun']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --norun", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option.

    # Other options
    #######################

#############################################################################