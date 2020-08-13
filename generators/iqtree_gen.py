#!/usr/bin/python
############################################################
# Generates commands for IQ-tree
############################################################

import sys, os, core, argparse

############################################################
# Options

parser = argparse.ArgumentParser(description="IQtree command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA alignments.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for tree runs.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to IQtree. Default: iqtree", default="iqtree");
parser.add_argument("-genetrees", dest="genetrees", help="Provide pre-estimated gene trees for the concordance analysis.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options

parser.add_argument("-part", dest="part", help="SLURM partition option.", default=False);
parser.add_argument("-tasks", dest="tasks", help="SLURM --ntasks option.", type=int, default=1);
parser.add_argument("-cpus", dest="cpus", help="SLURM --cpus-per-task option.", type=int, default=1);
parser.add_argument("-mem", dest="mem", help="SLURM --mem option.", type=int, default=0);
# SLURM options

args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory must be defined with -i.");
args.input = os.path.abspath(args.input);

if args.genetrees and not os.path.isfile(args.genetrees):
    sys.exit( " * Error 2: Cannot find gene tree file (-genetrees).");

if not args.name:
    name = core.getRandStr();
else:
    name = args.name;

if not args.output:
    sys.exit( " * Error 3: An output directory must be defined with -o.");

args.output = os.path.abspath(args.output);
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 4: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

if not args.part:
    sys.exit( " * Error 5: -part must be defined as a valid node partition on your clutser.");
if args.tasks < 1:
    sys.exit( " * Error 6: -tasks must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 7: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 8: -mem must be a positive integer.");
# SLURM option error checking

pad = 26
cwd = os.getcwd();
# Job vars

submit_file = os.path.join(cwd, "submit", name + ".sh");
# Job files

##########################
# Reporting run-time info for records.

core.runTime("# IQtree command generator");
core.PWS("# IO OPTIONS");
core.PWS(core.spacedOut("# Input directory:", pad) + args.input);
if args.genetrees:
    core.PWS(core.spacedOut("# Pre-estimated gene trees:", pad) + args.genetrees);
if not args.name:
    core.PWS("# -n not specified --> Generating random string for job name");
core.PWS(core.spacedOut("# Job name:", pad) + name);
core.PWS(core.spacedOut("# Output directory:", pad) + args.output);
if args.overwrite:
    core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.");
if not os.path.isdir(args.output):
    core.PWS("# Creating output directory.");
    os.system("mkdir " + args.output);
core.PWS("# ----------");
core.PWS("# SLURM OPTIONS");
core.PWS(core.spacedOut("# Submit file:", pad) + submit_file);
core.PWS(core.spacedOut("# SLURM partition:", pad) + args.part);
core.PWS(core.spacedOut("# SLURM ntasks:", pad) + str(args.tasks));
core.PWS(core.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus));
core.PWS(core.spacedOut("# SLURM mem:", pad) + str(args.mem));
core.PWS("# ----------");
core.PWS("# BEGIN CMDS");
    
##########################
# Generating the submit script.

if not args.genetrees:
    gene_tree_block = '''
echo "GENE TREES"
locidir="$treedir/loci/"
mkdir $locidir
cd $locidir
lociname="{name}-loci"
iqtree -S $alndir --prefix $lociname -T 1
'''
    gene_tree_block = gene_tree_block.format(name=name);

    locitrees = "$locidir/$lociname.treefile";

else:
    gene_tree_block = "";
    locitrees = args.genetrees;

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


alndir="{alndir}"
treedir="{treedir}"
mkdir $treedir

{gene_tree_block}

echo "CONCAT"
concatdir="$treedir/concat/"
mkdir $concatdir
cd $concatdir
concatname="{name}-concat"
iqtree -p $alndir --prefix $concatname -B 1000 -T {cpus}

echo "CONCORD"
cd $treedir
concorddir="$treedir/concord/"
mkdir $concorddir
cd $concorddir
locitrees="{locitrees}"
concattree="$concatdir/$concatname.treefile"
concordname="{name}-concord"
iqtree -t $concattree --gcf $locitrees -p $alndir --df-tree --scf 100 --cf-verbose --prefix $concordname -T 1
'''

    sfile.write(submit.format(name=name, partition=args.part, tasks=args.tasks, cpus=args.cpus, mem=args.mem, alndir=args.input, treedir=args.output, gene_tree_block=gene_tree_block, locitrees=locitrees));

##########################