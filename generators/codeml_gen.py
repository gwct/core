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

noisy = 3
verbose = 0
runmode = -2

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
aaRatefile = codeml/dat/wag.dat
model = 2

NSsites = 0

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
'''

############################################################
# Options

parser = argparse.ArgumentParser(description="codeml command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA alignments .", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to codeml. Default: codeml", default="codeml");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--outname", dest="outname", help="Use the end of the output directory path as the job name.", action="store_true", default=False);
# IO options

parser.add_argument("-tree", dest="tree", help="The species tree to use.", default=False);
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
if args.outname:
    name = os.path.basename(args.output);
else:
    args.output = args.output + "-" + name + "/";
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

if not args.tree or not os.path.isfile(args.tree):
    sys.exit(" * Error CO1: A -tree must be specified.");
args.tree = os.path.abspath(args.tree);
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

output_file = os.path.join(cwd, "jobs", "codeml_cmds_" + name + ".sh");
submit_file = os.path.join(cwd, "submit", "codeml_submit_" + name + ".sh");
#logdir = os.path.join(args.output, "logs");
# Job files

##########################
# Reporting run-time info for records.

with open(output_file, "w") as outfile:
    core.runTime("#!/bin/bash\n# codeml command generator", outfile);
    core.PWS("# IO OPTIONS", outfile);
    core.PWS(core.spacedOut("# Input directory:", pad) + args.input, outfile);
    if args.outname:
        core.PWS(core.spacedOut("# --outname:", pad) + "Using end of output directory path as job name.", outfile);
    if not args.name and not args.outname:
        core.PWS("# -n not specified --> Generating random string for job name", outfile);
    core.PWS(core.spacedOut("# Job name:", pad) + name, outfile);
    core.PWS(core.spacedOut("# Output directory:", pad) + args.output, outfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", outfile);
    if not os.path.isdir(args.output):
        core.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + args.output);
    # core.PWS(core.spacedOut("# Logfile directory:", pad) + logdir, outfile);
    # if not os.path.isdir(logdir):
    #     core.PWS("# Creating logfile directory.", outfile);
    #     os.system("mkdir " + logdir);
    core.PWS(core.spacedOut("# Job file:", pad) + output_file, outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# codeml OPTIONS", outfile);
    core.PWS(core.spacedOut("# Tree:", pad) + args.tree, outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# SLURM OPTIONS", outfile);
    core.PWS(core.spacedOut("# Submit file:", pad) + submit_file, outfile);
    core.PWS(core.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    core.PWS(core.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    core.PWS(core.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    core.PWS(core.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    core.PWS("# ----------", outfile);
    
##########################
# Generating the commands in the job file.

    num_skipped = 0;
    for f in os.listdir(args.input):
        base_input = os.path.splitext(f)[0];
        cur_infile = os.path.join(args.input, f);
        titles, seqs = core.fastaGetLists(cur_infile);
        seqs = [ s.lower() for s in seqs ];

        prem_stop = False
        for seq in seqs:
            if coreseq.premStopCheck(seq, allowlastcodon=True):
                prem_stop = True;
                break;

        if prem_stop:
            num_skipped += 1;
            continue;
        # if seqs.count(seqs[0]) == len(seqs):
        #     num_skipped += 1;
        #     continue;

        cur_outdir = os.path.join(args.output, base_input);
        if not os.path.isdir(cur_outdir):
            os.system("mkdir " + cur_outdir);

        cur_ctlfile = os.path.join(cur_outdir, "codeml.ctl");
        cur_outfile = os.path.join(cur_outdir, base_input + "-codeml-" + name + ".txt");
        cur_logfile = os.path.join(cur_outdir, base_input + "-codeml-" + name + ".log");

        with open(cur_ctlfile, "w") as ctlfile:
            ctlfile.write(ctlfile_template.format(infile=cur_infile, treefile=args.tree, outfile=cur_outfile));

        codeml_cmd = "cd " + cur_outdir + "; " + args.path + " codeml.ctl";

        outfile.write(codeml_cmd + "\n");

    core.PWS("# Num skipped because of premature stop codons: " + str(num_skipped));

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