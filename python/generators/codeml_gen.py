#!/usr/bin/python
############################################################
# Generates commands for codeml
############################################################

import sys, os, core, coreseq, argparse

############################################################
# Control file template

ctlfile_template = '''seqfile = {infile}
treefile = {treefile}
outfile = {outfile}

model = 0
NSsites = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
aaRatefile = codeml/dat/wag.dat

icode = 0
fix_kappa = 0
kappa = 3
fix_omega = 0
omega = 1

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
getSE = 0
Small_Diff = .5e-6

noisy = 3
verbose = 0
'''

############################################################
# Options

parser = argparse.ArgumentParser(description="codeml command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA alignments .", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to codeml. Default: codeml", default="codeml");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--stopcheck", dest="stopcheck", help="Set to check the input sequences for premature stop codons.", action="store_true", default=False);
# IO options

parser.add_argument("-tree", dest="tree", help="The species tree to use.", default=False);
parser.add_argument("-genetrees", dest="genetrees", help="A file containing a list of gene trees to use for each analysis, one per line. Must be in same order as alignment files.", default=False);
# Program options

parser.add_argument("-part", dest="part", help="SLURM partition option.", default=False);
parser.add_argument("-tasks", dest="tasks", help="SLURM --ntasks option.", type=int, default=1);
parser.add_argument("-cpus", dest="cpus", help="SLURM --cpus-per-task option.", type=int, default=1);
parser.add_argument("-mem", dest="mem", help="SLURM --mem option.", type=int, default=0);
# SLURM options

args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory must be defined with -i.");
args.input = os.path.abspath(args.input);

if not args.name:
    name = core.getRandStr();
else:
    name = args.name;

if not args.output:
    sys.exit( " * Error 2: An output directory must be defined with -o.");

args.output = os.path.abspath(args.output);
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

if args.tree and args.genetrees:
    sys.exit(" * Error CO1: Only one of -tree or -genetrees should be provided.");
if args.tree:
    if not os.path.isfile(args.tree):
        sys.exit(" * Error CO2: -tree: invalid file name.");
    treefile = os.path.abspath(args.tree);
elif args.genetrees:
    if not os.path.isfile(args.genetrees):
        sys.exit(" * Error CO3: -genetrees: invalid file name.");
    treefile = os.path.abspath(args.genetrees);
else:
    sys.exit(" * Error CO4: At least one of -tree or -genetrees must be provided.");
# codeml error checking

if not args.part:
    sys.exit( " * Error 4: -part must be defined as a valid node partition on your clutser.");
if args.tasks < 1:
    sys.exit( " * Error 5: -tasks must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 6: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 7: -mem must be a positive integer.");
# SLURM option error checking

pad = 26
cwd = os.getcwd();
# Job vars

output_file = os.path.join(cwd, "jobs", name + ".sh");
submit_file = os.path.join(cwd, "submit", name + ".sh");
#logdir = os.path.join(args.output, "logs");
# Job files

##########################
# Reporting run-time info for records.

with open(output_file, "w") as outfile:
    core.runTime("#!/bin/bash\n# codeml command generator", outfile);
    core.PWS("# IO OPTIONS", outfile);
    core.PWS(core.spacedOut("# Input directory:", pad) + args.input, outfile);
    if not args.name:
        core.PWS("# -n not specified --> Generating random string for job name", outfile);
    core.PWS(core.spacedOut("# Job name:", pad) + name, outfile);
    core.PWS(core.spacedOut("# Output directory:", pad) + args.output, outfile);
    if args.stopcheck:
        core.PWS(core.spacedOut("# --stopcheck set:", pad) + "Checking input sequences for premature stop codons.", outfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", outfile);
    if not os.path.isdir(args.output):
        core.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + args.output);
    core.PWS(core.spacedOut("# Job file:", pad) + output_file, outfile);
    core.PWS("# ----------------", outfile);

    core.PWS("# codeml OPTIONS", outfile);
    if args.tree:
        core.PWS(core.spacedOut("# Using single species tree:", pad) + treefile, outfile);
    elif args.genetrees:
        core.PWS(core.spacedOut("# Using gene trees:", pad) + treefile, outfile);
    core.PWS("# ----------------", outfile);

    core.PWS("# SLURM OPTIONS", outfile);
    core.PWS(core.spacedOut("# Submit file:", pad) + submit_file, outfile);
    core.PWS(core.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    core.PWS(core.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    core.PWS(core.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    core.PWS(core.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    core.PWS("# ----------------", outfile);
    core.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.

    aligns = sorted([ f for f in os.listdir(args.input) if f.endswith(".fa") ]);
    # Read and sort the alignment file names

    if args.genetrees:
        trees = list(filter(None, open(treefile, "r").read().split("\n")));
        assert len(aligns) == len(trees), "\n * Error 8: Unequal number of alignments (" + str(len(aligns)) + ") and gene trees (" + str(len(trees)) + ").";
    # Read all the gene trees if -genetrees is set
    else:
        trees = open(treefile, "r").read().split("\n");
    # Read the species tree if -tree is set

    num_skipped = 0;
    for i in range(len(aligns)):
        cur_align = aligns[i];
        if args.genetrees:
            cur_tree = trees[i];
        # Assign the current gene tree if -genetrees is set
        else:
            cur_tree = trees[0];
        # Simply use the species tree if -tree is set

        base_input = os.path.splitext(cur_align)[0];
        cur_infile = os.path.join(args.input, cur_align);
        # Get full paths to the input files

        seq_dict = core.fastaGetDict(cur_infile);
        if args.stopcheck:
            prem_stop_flag = False
            for title in seq_dict:
                prem_stop, new_seq = coreseq.premStopCheck(seq_dict[title], allowlastcodon=True, rmlast=True);
                if prem_stop:
                    prem_stop_flag = True;
                seq_dict[title] = new_seq;

            if prem_stop_flag:
                outfile.write(" # Premature stop found. Skipping: " + cur_align + "\n");
                num_skipped += 1;
                continue;
        # Check the alignment for premature stop codons (which PAML hangs on)

        cur_outdir = os.path.join(args.output, base_input);
        if not os.path.isdir(cur_outdir):
            os.system("mkdir " + cur_outdir);
        # Make the output directory for this alignment

        new_treefile = os.path.join(cur_outdir, "codeml.tre");
        with open(new_treefile, "w") as treefile:
            treefile.write(cur_tree);
        # Make the tree file for this alignment

        new_seqfile = os.path.join(cur_outdir, "codeml.fa");
        with open(new_seqfile, "w") as seqfile:
            for title in seq_dict:
                seqfile.write(title + "\n");
                seqfile.write(seq_dict[title] + "\n");
        # Write the sequences for this alignment

        cur_ctlfile = os.path.join(cur_outdir, "codeml.ctl");
        cur_outfile = os.path.join(cur_outdir, "codeml.out");
        #cur_logfile = os.path.join(cur_outdir, base_input + "-codeml.log");
        # Get the control and output file names

        with open(cur_ctlfile, "w") as ctlfile:
            ctlfile.write(ctlfile_template.format(infile=new_seqfile, treefile=new_treefile, outfile=cur_outfile));
         # Write the control file

        codeml_cmd = "cd " + cur_outdir + "; " + args.path + " codeml.ctl";
        outfile.write(codeml_cmd + "\n");
        # Construct and write the codeml command

    if args.stopcheck:
        core.PWS("# Num skipped because of premature stop codons: " + str(num_skipped), outfile);

##########################
# Generating the submit script.

with open(submit_file, "w") as sfile:
    submit = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --output={name}-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gregg.thomas@umontana.edu
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks={tasks}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=name, partition=args.part, tasks=args.tasks, cpus=args.cpus, mem=args.mem, output_file=output_file));

##########################