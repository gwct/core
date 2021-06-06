#!/usr/bin/python
############################################################
# Generates commands for the IQtree to make 1 gene tree at
# a time (since I can't get the -T option to work),
# parallelized with GNU Parallel.
# Also makes SLURM script for concatenation and concordance
# analysis.
############################################################

import sys, os, core, argparse

############################################################
# Options

parser = argparse.ArgumentParser(description="MACSE command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA alignment files.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-b", dest="bootstrap", help="The number of bootstrap replicates to perform. Default: 0", default="0");
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to MACSE. Default: iqtree", default="iqtree");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--outname", dest="outname", help="Use the end of the output directory path as the job name.", action="store_true", default=False);
# IO options

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

if not args.name:
    name = core.getRandStr();
else:
    name = args.name;

if not args.output:
    sys.exit( " * Error 2: An output directory must be defined with -o.");

args.output = os.path.abspath(args.output);
if args.outname:
    name = os.path.basename(args.output);
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");

if args.bootstrap < 1000:
    sys.exit(" * Error 4: Bootstrap (-b) must be at least 1000.");
# IO option error checking

if not args.part:
    sys.exit( " * Error 5: -part must be defined as a valid node partition on your clutser.");
# if args.nodes < 1:
#     sys.exit( " * Error 6: -nodes must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 7: -tasks must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 8: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 9: -mem must be a positive integer.");
# SLURM option error checking

pad = 26
cwd = os.getcwd();
# Job vars

locitree_dir = os.path.join(args.output, "loci");
locitree_file = os.path.join(args.output, "loci.treefile");
concat_dir = os.path.join(args.output, "concat");
concord_dir = os.path.join(args.output, "concord");
# Output sub-directories

loci_output_file = os.path.join(cwd, "jobs", name + "_loci.sh");
loci_submit_file = os.path.join(cwd, "submit", name + "_loci_submit.sh");
loci_logdir = os.path.join(args.output, "logs");

concat_submit_file = os.path.join(cwd, "submit", name + "_concat_submit.sh");
# Job files

##########################
# Reporting run-time info for records.

with open(loci_output_file, "w") as outfile:
    core.runTime("#!/bin/bash\n# IQtree command generator", outfile);
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
    core.PWS(core.spacedOut("# Loci tree directory:", pad) + locitree_dir, outfile);
    if not os.path.isdir(locitree_dir):
        core.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + locitree_dir);
    core.PWS(core.spacedOut("# Loci tree file:", pad) + locitree_file, outfile);
    core.PWS(core.spacedOut("# Concatenation directory:", pad) + concat_dir, outfile);
    if not os.path.isdir(concat_dir):
        core.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + concat_dir);
    core.PWS(core.spacedOut("# Concordance directory:", pad) + concord_dir, outfile);
    if not os.path.isdir(concord_dir):
        core.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + concord_dir);
    core.PWS(core.spacedOut("# Logfile directory:", pad) + loci_logdir, outfile);
    if not os.path.isdir(loci_logdir):
        core.PWS("# Creating logfile directory.", outfile);
        os.system("mkdir " + loci_logdir);
    core.PWS(core.spacedOut("# Job file:", pad) + loci_output_file, outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# SLURM OPTIONS", outfile);
    core.PWS(core.spacedOut("# Loci submit file:", pad) + loci_submit_file, outfile);
    core.PWS(core.spacedOut("# Concat submit file:", pad) + concat_submit_file, outfile);
    core.PWS(core.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    core.PWS(core.spacedOut("# SLURM nodes:", pad) + str(args.nodes), outfile);
    core.PWS(core.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    core.PWS(core.spacedOut("# SLURM ntasks-per-node:", pad) + str(args.tpn), outfile);
    core.PWS(core.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    core.PWS(core.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.
    skipped = 0;
    for f in sorted(os.listdir(args.input)):
        if not f.endswith(".fa") and not f.endswith(".phy"):
            continue;

        base_input = os.path.splitext(f)[0];
        cur_infile = os.path.join(args.input, f);

        cur_outdir = os.path.join(locitree_dir, base_input);
        if not os.path.isdir(cur_outdir):
            os.makedirs(cur_outdir);

        cur_out_prefix = os.path.join(cur_outdir, base_input);
        cur_logfile = os.path.join(loci_logdir, base_input + "-iqtree.log");

        iqtree_cmd = args.path + " -s " + cur_infile + " --prefix " + cur_out_prefix;
        if args.bootstrap > 0:
            iqtree_cmd += " -B " + str(args.bootstrap);
        iqtree_cmd += " -T 1 > " + cur_logfile + " 2>&1";

        outfile.write(iqtree_cmd + "\n");

    cat_cmd = "cat " + locitree_dir + "/*/*.treefile > " + locitree_file;
    outfile.write(cat_cmd + "\n");
    # Loci commands

    concat_log = os.path.join(concat_dir, "concat-terminal.log");
    concat_prefix = os.path.join(concat_dir, args.name);
    concat_cmd = "iqtree -p " + args.input + " --prefix " + concat_prefix + " -B 1000 -T " + str(args.cpus) + " &> " + concat_log;
    # Concatenation command

    concord_log = os.path.join(concord_dir, "concord-terminal.log");
    concat_file = os.path.join(concat_dir, args.name + ".treefile");
    concord_prefix = os.path.join(concord_dir, args.name);
    concord_cmd = "iqtree -t " + concat_file + " --gcf " + locitree_file + " -p " + args.input + " --scf 100 --cf-verbose --prefix " + concord_prefix + " -T 1";
    # Concordance command

    core.PWS("# ----------", outfile);
    core.PWS(core.spacedOut("# Files skipped: ", pad) + str(skipped), outfile);
    core.PWS("# Writing concat commands to " + concat_submit_file, outfile);
    core.PWS("# " + concat_cmd, outfile);
    core.PWS("# " + concord_cmd, outfile);

##########################
# Generating the loci submit script.

with open(loci_submit_file, "w") as sfile:
    submit = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --output={name}-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gregg.thomas@umontana.edu
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={tasks}
#SBATCH --ntasks-per-node={tpn}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=name, partition=args.part, nodes=args.nodes, tasks=args.tasks, tpn=args.tpn, cpus=args.cpus, mem=args.mem, output_file=loci_output_file));

##########################

##########################
# Generating the concat submit script.

with open(concat_submit_file, "w") as sfile:
    submit = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --output={name}-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gregg.thomas@umontana.edu
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={tasks}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

{concat_cmd}
{concord_cmd}'''

    sfile.write(submit.format(name=name, partition=args.part, nodes=args.nodes, tasks=1, cpus=args.tasks, mem=args.mem, concat_cmd=concat_cmd, concord_cmd=concord_cmd));

##########################