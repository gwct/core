#!/usr/bin/python
############################################################
# Generates commands for the MAFFT alignment program
############################################################

import sys, os, core, argparse

############################################################
# Options

parser = argparse.ArgumentParser(description="MAFFT command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA files.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to MAFFT. Default: mafft", default="mafft");
parser.add_argument("--accurate", dest="accurate", help="Set for more accurate (but slower) alignments", action="store_true", default=False);
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
print(args.input);
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
# if args.outname:
#     name = os.path.basename(args.output);
# else:
#     args.output = args.output + "-" + name + "/";
if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# IO option error checking

if not args.part:
    sys.exit( " * Error 4: -part must be defined as a valid node partition on your clutser.");
if args.tasks < 1:
    sys.exit( " * Error 5: -tasks must be a positive integer.");
if not args.tpn:
    args.tpn = args.tasks;
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
logdir = os.path.join(args.output, "logs");
# Job files

##########################
# Reporting run-time info for records.

with open(output_file, "w") as outfile:
    core.runTime("#!/bin/bash\n# MAFFT command generator", outfile);
    core.PWS("# IO OPTIONS", outfile);
    core.PWS(core.spacedOut("# Input directory:", pad) + args.input, outfile);
    if args.outname:
        core.PWS(core.spacedOut("# --outname:", pad) + "Using end of output directory path as job name.", outfile);
    if not args.name:
        core.PWS("# -n not specified --> Generating random string for job name", outfile);
    core.PWS(core.spacedOut("# Job name:", pad) + name, outfile);
    core.PWS(core.spacedOut("# Output directory:", pad) + args.output, outfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", outfile);
    if not os.path.isdir(args.output):
        core.PWS("# Creating output directory.", outfile);
        os.system("mkdir " + args.output);
    core.PWS(core.spacedOut("# Logfile directory:", pad) + logdir, outfile);
    if not os.path.isdir(logdir):
        core.PWS("# Creating logfile directory.", outfile);
        os.system("mkdir " + logdir);
    if args.accurate:
        core.PWS("# Running more --accurate alignments with --localpair --maxiterate 1000 --adjustdirection --op 3 --ep 0.123", outfile);
    core.PWS(core.spacedOut("# Job file:", pad) + output_file, outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# SLURM OPTIONS", outfile);
    core.PWS(core.spacedOut("# Submit file:", pad) + submit_file, outfile);
    core.PWS(core.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    core.PWS(core.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    core.PWS(core.spacedOut("# SLURM ntasks-per-node:", pad) + str(args.tpn), outfile);
    core.PWS(core.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    core.PWS(core.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.

    for f in os.listdir(args.input):
        if f.endswith(".fa") or f.endswith(".fa.filtered"):


            base_input = os.path.splitext(f)[0];
            cur_infile = os.path.join(args.input, f);
            cur_outfile = os.path.join(args.output, base_input + "-mafft.fa");
            cur_logfile = os.path.join(logdir, base_input + "-mafft.log");

            #mafft_cmd = args.path + " --preservecase " + cur_infile + " 2> " + cur_logfile + " 1> " + cur_outfile;
            if args.accurate:
                mafft_cmd = args.path + " --adjustdirection --preservecase --localpair --maxiterate 1000 --op 3 --ep 0.123 "  + cur_infile + " 2> " + cur_logfile + " 1> " + cur_outfile;
            else:
                mafft_cmd = args.path + " --adjustdirection --preservecase " + cur_infile + " 2> " + cur_logfile + " 1> " + cur_outfile;

            outfile.write(mafft_cmd + "\n");

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
#SBATCH --ntasks-per-node={tpn}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=name, partition=args.part, nodes=args.nodes, tasks=args.tasks, tpn=args.tpn, cpus=args.cpus, mem=args.mem, output_file=output_file));

##########################