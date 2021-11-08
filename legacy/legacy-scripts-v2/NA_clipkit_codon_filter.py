#!/usr/bin/python
############################################################
# For Penn genomes, 06.2020
# Takes a log file from a clipkit run on amino acid sequence
# and removes corresponding sites from codon alignment.
############################################################

import sys, os, core, coreseq, argparse

############################################################
# Options

parser = argparse.ArgumentParser(description="ClipKit codon filter");
parser.add_argument("-i", dest="aa_input", help="Directory of ClipKit filtered amino acid alignments and log files.", default=False);
parser.add_argument("-c", dest="cds_input", help="Directory of unfiltered CDS alignments.", default=False);
parser.add_argument("-o", dest="output", help="Desired output directory for filtered CDS alignments.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.aa_input or not os.path.isdir(args.aa_input):
    sys.exit( " * Error 1: An input directory with ClipKit filtered amino acid sequences must be defined with -i.");
args.aa_input = os.path.abspath(args.aa_input);

if not args.cds_input or not os.path.isdir(args.cds_input):
    sys.exit( " * Error 2: An input directory with unfiltered CDS sequences must be defined with -c.");
args.cds_input = os.path.abspath(args.cds_input);

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

pad = 26
cwd = os.getcwd();
# Job vars

log_file = os.path.join("logs", name + ".log");
# Job files

##########################
# Reporting run-time info for records.

with open(log_file, "w") as logfile:
    core.runTime("# ClipKit CDS filter", logfile);
    core.PWS("# IO OPTIONS", logfile);
    core.PWS(core.spacedOut("# Input AA directory:", pad) + args.aa_input, logfile);
    core.PWS(core.spacedOut("# Input CDS directory:", pad) + args.cds_input, logfile);
    if not args.name:
        core.PWS("# -n not specified --> Generating random string for job name", logfile);
    core.PWS(core.spacedOut("# Job name:", pad) + name, logfile);
    core.PWS(core.spacedOut("# Output directory:", pad) + args.output, logfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", logfile);
    if not os.path.isdir(args.output):
        core.PWS("# Creating output directory.", logfile);
        os.system("mkdir " + args.output);
    core.PWS(core.spacedOut("# Log file:", pad) + log_file, logfile);
    core.PWS("# ----------------", logfile);
    
##########################
# Filtering CDS aligns

    core.PWS("# " + core.getDateTime() + " Beginning filter...", logfile);
    headers = "Align\tInitial length\tSites removed\tPercent sites removed\tPercent high\tFinal length\tShort aln\tStop codons";
    logfile.write(headers + "\n");
    skipped, seq_prem_stop, aln_prem_stop, num_high, num_short = 0,0,0,0,0;
    for f in os.listdir(args.aa_input):
        if not f.endswith(".fa.log"):
            continue;

        base_input = f.replace(".clipkit.fa.log", "");
        clip_log = os.path.join(args.aa_input, f);
        cds_aln = os.path.join(args.cds_input, base_input.replace("AA", "NT") + ".fa");

        cur_outfile = os.path.join(args.output, f.replace("AA", "NT").replace(".fa.log", ".fa"));
        # Get the current in and output files

        if not os.path.isfile(clip_log):
            core.PWS(cds_aln + "\t" + " Skipping -- could not find ClipKit log file:\t" + clip_log);
            skipped += 1;
        if not os.path.isfile(cds_aln):
            core.PWS(cds_aln + "\t" + " Skipping -- could not find CDS align file.", logfile);
            skipped += 1;
        # Make sure all files exist
        # print(f);
        # print(base_input);
        # print(clip_log);
        # print(cds_aln);
        # print(cur_outfile);
        # sys.exit();

        titles, seqs = core.fastaGetLists(cds_aln);
        # Read the CDS sequences

        codon_len = len(seqs[0]) / 3;
        # Get the total number of codons in the alignment

        rm_codons = [];
        for line in open(clip_log):
            line = line.strip().split(" ");
            if line[1] == "trim":
                rm_codons.append(int(line[0])-1);
        # Read the clipkit log to get the codons to be removed

        perc_rm = round(float((len(rm_codons)) / float(codon_len)) * 100.0, 2);
        # Count percent of sites removed

        codon_seqs = [];
        for seq in seqs:
            codon_list = [ seq[i:i+3] for i in range(0, len(seq), 3) ];
            codon_seqs.append(codon_list);
        # Convert each CDS sequence into a list of codons

        for site in sorted(rm_codons, reverse=True):
            for codon_seq in codon_seqs:
                del codon_seq[site];
        # Remove each site from each list of codons. Loop through sites in reverse to preserve indices.

        new_seqs = [];
        for codon_seq in codon_seqs:
            new_seqs.append("".join(codon_seq));
        # Convert from codon lists to sequence strings

        prem_stop_flag = False;
        for s in range(len(new_seqs)):
            stop, new_seqs[s] = coreseq.premStopCheck(new_seqs[s], allowlastcodon=True, rmlast=True);
            if stop:
                seq_prem_stop += 1;
                prem_stop_flag = True;
        if prem_stop_flag:
            aln_prem_stop += 1;
        # Check for premature stop codons and remove last codon if it is a stop

        filtered_codon_len = len(new_seqs[0]) / 3;
        # Get the total number of codons in the filtered alignment

        with open(cur_outfile, "w") as outfile:
            for i in range(len(titles)):
                outfile.write(titles[i] + "\n");
                outfile.write(new_seqs[i].replace("!", "N") + "\n");
        # Write filtered sequences to output file

        perc_high = "FALSE";
        if perc_rm > 20:
            num_high += 1;
            perc_high = "TRUE";

        short_aln = "FALSE";
        if filtered_codon_len < 50:
            num_short += 1;
            short_aln = "TRUE";

        outline = [cds_aln, str(int(codon_len)), str(len(rm_codons)), str(perc_rm), perc_high, str(int(filtered_codon_len)), short_aln, str(seq_prem_stop)];
        core.PWS("\t".join(outline), logfile);

    core.PWS("# ----------------", logfile);
    core.PWS(core.spacedOut("# Aligns with >20% codons removed:", 55) + str(num_high), logfile); 
    core.PWS(core.spacedOut("# Aligns shorter than 50 codons after filter", 55) + str(num_short), logfile); 
    core.PWS(core.spacedOut("# Aligns with premature stop:", 55) + str(aln_prem_stop), logfile);               
    core.PWS("# ----------------", logfile);


# java -jar ~/bin/macse_v2.03.jar -prog trimNonHomologousFragments -seq ENSMUSG00000042419-rodent-26-NEW.fa -out_NT -out_AA
# java -jar ~/bin/macse_v2.03.jar -prog alignSequences -seq ENSMUSG00000042419-rodent-26-NEW_NT.fa -out_NT -out_AA
# clipkit ENSMUSG00000042419-rodent-26-NEW_NT_AA.fa -m kpic-gappy -l
# this script
# java -jar ~/bin/macse_v2.03.jar -prog trimAlignment -seq ENSMUSG00000042419-rodent-26-NEW_NT_NT.fa.clipkit

