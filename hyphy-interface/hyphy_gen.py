#!/usr/bin/python
############################################################
# Generates commands and shell scripts for various PAML models
############################################################

import sys, os, re, argparse, lib.hpcore as hpcore, lib.hptree as hptree

############################################################
# Options

parser = argparse.ArgumentParser(description="codeml command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA alignments. Note: for -model anc-recon this should be a directory of Hyphy .json files.", default=False);
parser.add_argument("-m", dest="model", help="The model to run. Options: mg94, mg94-local, rm-dup, fel, busted, fubar, absrel, anc-recon, slac, relax. Default: mg94", default="mg94");
parser.add_argument("-s", dest="sep", help="The character to split the alignment filename on to obtain the gene ID. Default: none, use whole file name.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to codeml. Default: codeml", default="codeml");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--outname", dest="outname", help="Use the end of the output directory path as the job name.", action="store_true", default=False);
# IO options

parser.add_argument("-tree", dest="tree", help="The species tree to use.", default=False);
parser.add_argument("-genetrees", dest="genetrees", help="A directory containing gene trees for each locus (from iqtree_gt_gen.py).", default=False);
parser.add_argument("-targetclades", dest="target_clades", help="A file or directory of files containing subtrees of target clades.", default=False);
#parser.add_argument("-target", dest="target", help="A single or pair of species. The MRCA of the given species will be determined as the target lineage. Pairs should be separated by semi-colons and multiple targets separated by commas, e.g. 'targ1s1;targ1s2,targ2s1;targ2s2", default=False);
parser.add_argument("-tb", dest="testbranches", help="A comma delimited list of species or branches that make up the test branches for RELAX.", default=False);
parser.add_argument("-rb", dest="refbranches", help="A comma delimited list of species or branches that make up the reference branches for RELAX.", default=False);
# Program options

parser.add_argument("-part", dest="part", help="SLURM partition option.", default=False);
parser.add_argument("-nodes", dest="nodes", help="SLURM --nodes option.", default="1");
parser.add_argument("-tasks", dest="tasks", help="SLURM --ntasks option.", type=int, default=1);
parser.add_argument("-tpn", dest="tpn", help="SLURM --ntasks-per-node-option option.", default=False);
parser.add_argument("-cpus", dest="cpus", help="SLURM --cpus-per-task option.", type=int, default=1);
parser.add_argument("-mem", dest="mem", help="SLURM --mem option.", type=int, default=0);
# SLURM options

args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory must be defined with -i.");
args.input = os.path.abspath(args.input);

if args.model not in ["mg94", "mg94-local", "rm-dup", "fel", "busted", "fubar", "absrel", "anc-recon", "slac", "relax"]:
    sys.exit(" * Error 2: Model (-m) must be one of: mg94, mg94-local, rm-dup, fel, busted, fubar, absrel, slac, relax");

if args.target_clades:
    targets = {};
    if os.path.isfile(args.target_clades):
        targets.append(hptree.readTips(args.target_clades));
    elif os.path.isdir(args.target_clades):
        for target_file in os.listdir(args.target_clades):
            targets[target_file] = hptree.readTips(os.path.join(args.target_clades, target_file));
    else:
        sys.exit(" * Error 3: Target file/directory (-targetclades) not found!");
else:
    targets = False;
# Parse the targets.

if args.testbranches:
    tests = args.testbranches.replace(", ", ",").split(",");
    tests = list(set(tests));
    tests.sort();
else:
    tests = False;
# Parse the test branches.

if args.refbranches:
    refs = args.refbranches.replace(", ", ",").split(",");
    refs = list(set(refs));
    refs.sort();
else:
    refs = False;
# Parse the reference branches.

if not args.name:
    name = hpcore.getRandStr();
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

if args.model != "anc-recon":
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

if args.model == "relax":
    if not tests:
        sys.exit(" * Error 10: If running RELAX, -tb must be provided.");
    if not refs:
        sys.exit(" * Error 11: If running RELAX, -rb must be provided.");
# error checking

if not args.part:
    sys.exit( " * Error 12: -part must be defined as a valid node partition on your clutser.");
if args.tasks < 1:
    sys.exit( " * Error 13: -tasks must be a positive integer.");
if not args.tpn:
    args.tpn = args.tasks;
if args.cpus < 1:
    sys.exit( " * Error 14: -cpus must be a positive integer.");
if args.mem < 0:
    sys.exit( " * Error 15: -mem must be an integer >= 0.");
# SLURM option error checking

pad = 26
cwd = os.getcwd();
# Job vars

output_file = os.path.join(cwd, "jobs", name + ".sh");
submit_file = os.path.join(cwd, "submit", name + ".sh");
logdir = os.path.join(args.output, "logs");
treedir = os.path.join(args.output, "trees");
# Job files

##########################
# Reporting run-time info for records.

with open(output_file, "w") as outfile:
    hpcore.runTime("#!/bin/bash\n# HyPhy command generator", outfile);
    hpcore.PWS("# IO OPTIONS", outfile);
    hpcore.PWS(hpcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    hpcore.PWS(hpcore.spacedOut("# HyPhy model:", pad) + args.model, outfile);
    if args.outname:
        hpcore.PWS(hpcore.spacedOut("# --outname:", pad) + "Using end of output directory path as job name.", outfile);
    if not args.name and not args.outname:
        hpcore.PWS("# -n not specified --> Generating random string for job name", outfile);
    hpcore.PWS(hpcore.spacedOut("# Job name:", pad) + name, outfile);
    hpcore.PWS(hpcore.spacedOut("# Output directory:", pad) + args.output, outfile);
    hpcore.PWS(hpcore.spacedOut("# Log directory:", pad) + logdir, outfile);
    if args.overwrite:
        hpcore.PWS(hpcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", outfile);
    if not os.path.isdir(args.output):
        hpcore.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + args.output);
    if not os.path.isdir(logdir):
        hpcore.PWS("# Creating log directory.", outfile);
        os.system("mkdir " + logdir);
    if not os.path.isdir(treedir):
        hpcore.PWS("# Creating labeled tree directory.", outfile);
        os.system("mkdir " + treedir);
    hpcore.PWS(hpcore.spacedOut("# Job file:", pad) + output_file, outfile);
    hpcore.PWS("# ----------------", outfile);

    hpcore.PWS("# HyPhy OPTIONS", outfile);
    if args.tree:
        hpcore.PWS(hpcore.spacedOut("# Using single species tree:", pad) + tree_input, outfile);
    elif args.genetrees:
        hpcore.PWS(hpcore.spacedOut("# Using gene trees:", pad) + tree_input, outfile);
    if args.target_clades:
        hpcore.PWS(hpcore.spacedOut("# Target clade file/dir:", pad) + args.target_clades, outfile);
        hpcore.PWS(hpcore.spacedOut("# Num of target branches:", pad) + str(len(targets)), outfile);
    if args.testbranches:
        hpcore.PWS(hpcore.spacedOut("# Test branches:", pad) + ",".join(tests), outfile);
    if args.refbranches:
        hpcore.PWS(hpcore.spacedOut("# Reference branches:", pad) + ",".join(refs), outfile);
        if not os.path.isdir(treedir):
            hpcore.PWS("# Creating tree directory.", outfile);
            os.system("mkdir " + treedir);
    hpcore.PWS("# ----------------", outfile);

    hpcore.PWS("# SLURM OPTIONS", outfile);
    hpcore.PWS(hpcore.spacedOut("# Submit file:", pad) + submit_file, outfile);
    hpcore.PWS(hpcore.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    hpcore.PWS(hpcore.spacedOut("# SLURM nodes:", pad) + str(args.nodes), outfile);
    hpcore.PWS(hpcore.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    hpcore.PWS(hpcore.spacedOut("# SLURM ntasks-per-node:", pad) + str(args.tpn), outfile);
    hpcore.PWS(hpcore.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    hpcore.PWS(hpcore.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.

    file_dir = os.path.dirname(os.path.abspath(__file__));
    if args.model == "mg94":
        import lib.mg94 as mg94; 
        model_file = os.path.join(file_dir, "hyphy-analyses/FitMG94/FitMG94.bf");
        mg94.generate(args.input, tree_input, model_file, args.genetrees, args.sep, args.path, args.output, logdir, outfile);
    if args.model == "mg94-local":   
        import lib.mg94local as mg94local;
        model_file = os.path.join(file_dir, "hyphy-analyses/FitMG94/FitMG94.bf");
        mg94local.generate(args.input, tree_input, model_file, args.genetrees, args.path, args.output, logdir, outfile);
    if args.model == "rm-dup":   
        import lib.rmdup as rmdup;
        model_file = os.path.join(file_dir, "hyphy-analyses/remove-duplicates/remove-duplicates.bf");
        rmdup.generate(args.input, tree_input, model_file, args.genetrees, args.sep, args.path, args.output, logdir, outfile);

    if args.model == "fel":
        import lib.fel as fel;
        fel.generate(args.input, tree_input, args.genetrees, args.sep, args.path, args.output, logdir, outfile);
    if args.model == "busted":
        import lib.busted as busted;
        busted.generate(args.input, tree_input, args.genetrees, args.sep, args.path, args.output, logdir, outfile);
    if args.model == "fubar":
        import lib.fubar as fubar;
        fubar.generate(args.input, tree_input, args.genetrees, args.sep, args.path, args.output, logdir, outfile);
    if args.model == "absrel":
        import lib.absrel as absrel;
        absrel.generate(args.input, tree_input, args.genetrees, args.sep, targets, args.path, args.output, treedir, logdir, outfile);
    if args.model == "anc-recon":
        import lib.ancrecon as ancrecon;
        model_file = os.path.join(file_dir, "hyphy-analyses/AncestralSequences/AncestralSequences.bf");
        ancrecon.generate(args.input, model_file, args.path, args.output, logdir, outfile);
    if args.model == "slac":
        import lib.slac as slac;
        slac.generate(args.input, tree_input, args.genetrees, args.path, args.output, logdir, outfile);
    if args.model == "relax":
        import lib.relax as relax;
        relax.generate(args.input, tree_input, tests, refs, args.genetrees, args.path, args.output, logdir, outfile)


##########################
# Generating the submit script.

with open(submit_file, "w") as sfile:
    submit = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --output={name}-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gregg.thomas@umontana.edu
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={tasks}
#SBATCH --tasks-per-node={tpn}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

conda activate hyphyenv

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=name, partition=args.part, nodes=args.nodes, tasks=args.tasks, tpn=args.tpn, cpus=args.cpus, mem=args.mem, output_file=output_file));

##########################
