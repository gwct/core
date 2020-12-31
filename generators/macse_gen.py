#!/usr/bin/python
############################################################
# Generates commands for the MACSE alignment program
############################################################

import sys, os, core, argparse

############################################################
# Options

parser = argparse.ArgumentParser(description="MACSE command generator");
parser.add_argument("-i", dest="input", help="Directory of input FASTA files.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for aligned files. Job name (-n) will be appended to output directory name.", default=False);
parser.add_argument("-e", dest="expected", help="The expected number of species in each alignment file. Check for one-to-one alignments to only align sequences that retained all species after trimming.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("-p", dest="path", help="The path to MACSE. Default: java -jar ~/bin/macse_v2.03.jar", default="java -jar ~/bin/macse_v2.03.jar");
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
parser.add_argument("--outname", dest="outname", help="Use the end of the output directory path as the job name.", action="store_true", default=False);
# IO options

parser.add_argument("-part", dest="part", help="SLURM partition option.", default=False);
parser.add_argument("-nodes", dest="nodes", help="SLURM --nodes option.", type=int, default=1);
parser.add_argument("-tasks", dest="tasks", help="SLURM --ntasks option.", type=int, default=1);
parser.add_argument("-taskspernode", dest="taskspernode", help="SLURM --ntasks-per-node option.", type=int, default=1);
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

expected = False;
if args.expected:
    eflag = False;
    try:
        expected = int(args.expected);
    except:
        eflag = True;
    if eflag or expected < 1:
        sys.exit(" * Error 4: -e must be a positive integer.");
# IO option error checking

if not args.part:
    sys.exit( " * Error 5: -part must be defined as a valid node partition on your clutser.");
if args.nodes < 1:
    sys.exit( " * Error 6: -nodes must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 7: -tasks must be a positive integer.");
if args.taskspernode < 1:
    sys.exit( " * Error 8: -taskspernode must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 9: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 10: -mem must be a positive integer.");
# SLURM option error checking

pad = 26
cwd = os.getcwd();
# Job vars

output_file = os.path.join(cwd, "jobs", name + ".sh");
submit_file = os.path.join(cwd, "submit", name + ".sh");
aadir = os.path.join(args.output, "aa");
ntdir = os.path.join(args.output, "nt");
logdir = os.path.join(args.output, "logs");
# Job files

##########################
# Reporting run-time info for records.

with open(output_file, "w") as outfile:
    core.runTime("#!/bin/bash\n# MACSE command generator", outfile);
    core.PWS("# IO OPTIONS", outfile);
    core.PWS(core.spacedOut("# Input directory:", pad) + args.input, outfile);
    if expected:
        core.PWS(core.spacedOut("# Expected # of species:", pad) + str(expected), outfile);
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
    core.PWS(core.spacedOut("# Nucleotide directory:", pad) + ntdir, outfile);
    if not os.path.isdir(ntdir):
        core.PWS("# Creating nucleotide directory.", outfile);
        os.system("mkdir " + ntdir);
    core.PWS(core.spacedOut("# Amino acid directory:", pad) + aadir, outfile);
    if not os.path.isdir(aadir):
        core.PWS("# Creating amino acid directory.", outfile);
        os.system("mkdir " + aadir);
    core.PWS(core.spacedOut("# Logfile directory:", pad) + logdir, outfile);
    if not os.path.isdir(logdir):
        core.PWS("# Creating logfile directory.", outfile);
        os.system("mkdir " + logdir);
    core.PWS(core.spacedOut("# Job file:", pad) + output_file, outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# SLURM OPTIONS", outfile);
    core.PWS(core.spacedOut("# Submit file:", pad) + submit_file, outfile);
    core.PWS(core.spacedOut("# SLURM partition:", pad) + args.part, outfile);
    core.PWS(core.spacedOut("# SLURM nodes:", pad) + str(args.nodes), outfile);
    core.PWS(core.spacedOut("# SLURM ntasks:", pad) + str(args.tasks), outfile);
    core.PWS(core.spacedOut("# SLURM cpus-per-task:", pad) + str(args.cpus), outfile);
    core.PWS(core.spacedOut("# SLURM mem:", pad) + str(args.mem), outfile);
    core.PWS("# ----------", outfile);
    core.PWS("# BEGIN CMDS", outfile);
    
##########################
# Generating the commands in the job file.
    skipped = 0;
    for f in os.listdir(args.input):
        if not f.endswith(".fa"):
            continue;

        base_input = os.path.splitext(f)[0];
        cur_infile = os.path.join(args.input, f);

        if expected:
            inseqs = core.fastaGetDict(cur_infile);
            if len(inseqs) != expected:
                outfile.write("# Expected number of species not found... skipping: " + cur_infile + "\n");
                skipped += 1;
                continue;


        cur_outfile_nt = os.path.join(ntdir, base_input + ".NT.macse.fa");
        cut_outfile_aa = os.path.join(aadir, base_input + ".AA.macse.fa");
        cur_logfile = os.path.join(logdir, base_input + "-macse.log");

        macse_cmd = args.path + " -prog alignSequences -seq '" + cur_infile + "' -out_NT '" + cur_outfile_nt + "' -out_AA '" + cut_outfile_aa + "' > " + cur_logfile + " 2>&1";
        
        #macse_cmd = args.path + " -prog alignSequences -seq '" + cur_infile + "' -stop 200 -out_NT '" + cur_outfile_nt + "' -out_AA '" + cut_outfile_aa + "' > " + cur_logfile + " 2>&1";
        # Do not allow stop codons

        outfile.write(macse_cmd + "\n");

    core.PWS("# ----------", outfile);
    core.PWS(core.spacedOut("# Files skipped: ", pad) + str(skipped), outfile);

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
#SBATCH --ntasks-per-node={ntaskspernode}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=name, partition=args.part, nodes=args.nodes, tasks=args.tasks, ntaskspernode=args.taskspernode, cpus=args.cpus, mem=args.mem, output_file=output_file));

##########################