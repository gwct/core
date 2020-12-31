#!/usr/bin/python
############################################################
# Generates commands and shell scripts for various PAML models
############################################################

import sys, os, re, argparse, lib.pamlcore as pcore

############################################################
# Options

parser = argparse.ArgumentParser(description="codeml command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA alignments .", default=False);
parser.add_argument("-m", dest="model", help="The PAML model that was used to generate the files in -i. Options: m1, cmc. Default: m1", default="m1");
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to codeml. Default: codeml", default="codeml");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--outname", dest="outname", help="Use the end of the output directory path as the job name.", action="store_true", default=False);
# IO options

parser.add_argument("-tree", dest="tree", help="The species tree to use.", default=False);
parser.add_argument("-genetrees", dest="genetrees", help="A directory containing gene trees for each locus (from iqtree_gt_gen.py).", default=False);
parser.add_argument("-targetclade", dest="target_clade", help="For clade model C (-m cmc) and Model 2 (-m m2): A comma delimited list of species that make up the clade descending from the target branch. Only alignments and gene trees with the target clade will be used.", default=False);
parser.add_argument("--anc", dest="anc_recon", help="Set to also have codeml perform ancestral sequence reconstruction.", action="store_true", default=False);
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

if args.model not in ["m1", "m2", "cmc"]:
    sys.exit(" * Error 2: Model (-m) must be one of: m1, m2, cmc");

if args.model in ["m2", "cmc"]:
    if not args.target_clade:
        sys.exit(" * Error 3: With models m2 and cmc a -targetclade must be given.");

    targets = args.target_clade.replace(", ", ",").split(",");
    targets = set(targets);
# Parse the targets for CMC.

if not args.name:
    name = pcore.getRandStr();
else:
    name = args.name;

if not args.output:
    sys.exit( " * Error 4: An output directory must be defined with -o.");

args.output = os.path.abspath(args.output);
if args.outname:
    name = os.path.basename(args.output);
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 5: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

if args.tree and args.genetrees:
    sys.exit(" * Error 6: Only one of -tree or -genetrees should be provided.");
if args.tree:
    if not os.path.isfile(args.tree):
        sys.exit(" * Error 7: -tree: invalid file name.");
    tree_input = os.path.abspath(args.tree);
elif args.genetrees:
    if not os.path.isdir(args.genetrees):
        sys.exit(" * Error 8: -genetrees: invalid directory name.");
    tree_input = os.path.abspath(args.genetrees);
else:
    sys.exit(" * Error 9: At least one of -tree or -genetrees must be provided.");
anc_recon_setting = "0";
if args.anc_recon:
    anc_recon_setting = "1";
# codeml error checking

if not args.part:
    sys.exit( " * Error 10: -part must be defined as a valid node partition on your clutser.");
if args.tasks < 1:
    sys.exit( " * Error 11: -tasks must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 12: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 13: -mem must be a positive integer.");
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
    pcore.runTime("#!/bin/bash\n# codeml command generator", outfile);
    pcore.PWS("# IO OPTIONS", outfile);
    pcore.PWS(pcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    pcore.PWS(pcore.spacedOut("# PAML model:", pad) + args.model, outfile);
    if args.outname:
        pcore.PWS(pcore.spacedOut("# --outname:", pad) + "Using end of output directory path as job name.", outfile);
    if not args.name and not args.outname:
        pcore.PWS("# -n not specified --> Generating random string for job name", outfile);
    pcore.PWS(pcore.spacedOut("# Job name:", pad) + name, outfile);
    pcore.PWS(pcore.spacedOut("# Output directory:", pad) + args.output, outfile);
    if args.overwrite:
        pcore.PWS(pcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", outfile);
    if not os.path.isdir(args.output):
        pcore.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + args.output);
    pcore.PWS(pcore.spacedOut("# Job file:", pad) + output_file, outfile);
    pcore.PWS("# ----------------", outfile);

    pcore.PWS("# codeml OPTIONS", outfile);
    if args.tree:
        pcore.PWS(pcore.spacedOut("# Using single species tree:", pad) + tree_input, outfile);
    elif args.genetrees:
        pcore.PWS(pcore.spacedOut("# Using gene trees:", pad) + tree_input, outfile);
    if args.model in ["m2", "cmc"]:
        pcore.PWS(pcore.spacedOut("# Target clade:", pad) + ",".join(targets), outfile);
    if args.anc_recon:
        pcore.PWS(pcore.spacedOut("# Performing ancestral sequence reconstruction.", pad), outfile);
    pcore.PWS("# ----------------", outfile);

    pcore.PWS("# SLURM OPTIONS", outfile);
    pcore.PWS(pcore.spacedOut("# Submit file:", pad) + submit_file, outfile);
    pcore.PWS(pcore.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    pcore.PWS(pcore.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    pcore.PWS(pcore.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    pcore.PWS(pcore.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    pcore.PWS("# ----------------", outfile);
    pcore.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.

    if args.model == "m1":
        import lib.m1 as m1;
        m1.generate(args.input, tree_input, args.genetrees, anc_recon_setting, args.path, args.output, outfile);
    if args.model == "cmc":
        import lib.cmc as cmc;
        cmc.generate(args.input, tree_input, args.genetrees, anc_recon_setting, targets, args.path, args.output, outfile);
    if args.model == "m2":
        import lib.m2 as m2;
        m2.generate(args.input, tree_input, args.genetrees, anc_recon_setting, targets, args.path, args.output, outfile);


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