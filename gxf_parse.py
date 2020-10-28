#!/home/gt156213e/anaconda3/envs/main/bin/python
########################################################################################
# This script converts GFF or GTF files to a more sensible tab delimited file.
#
# Dependencies: core
#
# Gregg Thomas, Fall 2020
########################################################################################

import sys, os, re, argparse
from collections import defaultdict
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

########################################################################################
# Functions


########################################################################################
# Main

parser = argparse.ArgumentParser(description="GTF/GFF to generic tab conversion.");
parser.add_argument("-i", dest="input", help="The input GTF or GFF formatted file.", default=False);
parser.add_argument("-s", dest="source", help="The source of the file. Currently supported sources: 'ensembl' and 'maker'", default=False);
parser.add_argument("-p", dest="prefix", help="The feature ID prefix to retrieve gene/transcript/exon IDs.", default=False);
parser.add_argument("-o", dest="output", help="The output tab delmited file.", default=False);
parser.add_argument("--nogenes", dest="no_genes", help="Set to count the number of reads in each file.", default=False, action="store_true");
parser.add_argument("--notranscripts", dest="no_transcripts", help="Set to calculate the average read length in each file.", default=False, action="store_true");
parser.add_argument("--noexons", dest="no_exons", help="Set to count the base compositions in each file.", default=False, action="store_true");
parser.add_argument("--overwrite", dest="overwrite", help="Set this to indicate you wish to overwrite files specified by -outcsv and -outtxt if they already exist. WARNING: This means the original contents of the file will be deleted.", default=False, action="store_true");
#parser.add_argument("--header", dest="header", help="Set to extract the header info for each file.", default=False, action="store_true");
args = parser.parse_args();
# Input option definitions.

if not args.input or not os.path.isfile(args.input):
    sys.exit(core.errorOut(1, "Cannot find input file (-i)."));

if not args.source:
    sys.exit(core.errorOut(2, "Please specify either ensembl or maker as the source of the input file (-s)."));

if not args.prefix:
    sys.exit(core.errorOut(3, "Please specify a feature ID prefix (i.e. ENSMUS) for the current file (-p)."));

if not args.output:
    sys.exit(core.errorOut(4, "Please specify the name of an output file (-o)."));
elif os.path.isfile(args.output) and not args.overwrite:
   sys.exit(core.errorOut(5, "Output file (-o) already exists! Explicity specify --overwrite to overwrite it."));
# I/O parsing and error checking.

if args.source not in ["maker", "ensembl"]:
    sys.exit(core.errorOut(6, "Source (-s) of file must be: 'maker' or 'ensembl'"));

if args.source == "maker":
    transcript_str = "mRNA";
elif args.source == "ensembl":
    transcript_str = "transcript";

pad = 25;
with open(args.output, "w") as outfile:
    core.runTime("# GFF/GTF parsing and conversion.", outfile);
    core.PWS(core.spacedOut("# Input file:", pad) + args.input, outfile);
    core.PWS(core.spacedOut("# Input source:", pad) + args.source, outfile);
    core.PWS(core.spacedOut("# Feature ID prefix:", pad) + args.prefix, outfile);
    core.PWS(core.spacedOut("# Output file:", pad) + args.output, outfile);
    if args.overwrite:
        core.PWS(core.spacedOut("# --overwrite set:", pad) + "Overwriting previous output file.", outfile);
    if args.no_genes:
        core.PWS(core.spacedOut("# --nogenes set:", pad) + "Not writing gene coordinates and IDs.", outfile);
    if args.no_transcripts:
        core.PWS(core.spacedOut("# --notranscripts set:", pad) + "Not writing transcript coordinates and IDs.", outfile);
    if args.no_exons:
        core.PWS(core.spacedOut("# --noexons set:", pad) + "Not writing exon coordinates and IDs.", outfile);
    core.PWS("# ----------------", outfile);

    core.PWS("# " + core.getDateTime() + " Detecting compression...", outfile);
    reader = core.getFileReader(args.input);
    if reader == open:
        line_reader = core.readLine;
        read_mode = "r";
    if reader != open:
        line_reader = core.readGzipLine
        read_mode = "rb";

    core.PWS("# " + core.getDateTime() + " Reading genes...", outfile);
    genes = {};
    for line in reader(args.input, read_mode):
        #print(line);
        line = line_reader(line);
        if "##FASTA" in line[0]:
            break;
        if line[0][0] == "#":
            continue;
        feature_type, chrome, start, end, strand, feature_info = line[2], line[0], line[3], line[4], line[6], line[8];

        if feature_type == "gene":
            if args.source == "maker":
                gid = re.findall(args.prefix + '_G[\d]+', feature_info)[0];
            elif args.source == "ensembl":
                gid = re.findall(args.prefix + 'G[\d]+', feature_info)[0];
            genes[gid] = { 'chr' : chrome, 'start' : start, 'end' : end, 'strand' : strand, 'transcripts' : {} };
    core.PWS(core.spacedOut("# Genes read:", pad) + str(len(genes)), outfile);
    core.PWS("# ----------------", outfile);
    # Read genes.

    if not args.no_transcripts or not args.no_exons:
        core.PWS("# " + core.getDateTime() + " Reading transcripts...", outfile);

        num_transcripts = 0;
        for line in reader(args.input, read_mode):
            #print(line);
            line = line_reader(line);
            if "##FASTA" in line[0]:
                break;
            if line[0][0] == "#":
                continue;
            feature_type, chrome, start, end, strand, feature_info = line[2], line[0], line[3], line[4], line[6], line[8];

            if feature_type == transcript_str:
                if args.source == "maker":
                    gid = re.findall(args.prefix + '_G[\d]+', feature_info)[0];
                    tid = re.findall(args.prefix + '_T[\d]+_mRNA-[\d]+', feature_info)[0];
                    transcript_type = "";
                elif args.source == "ensembl":
                    gid = re.findall(args.prefix + 'G[\d]+', feature_info)[0];
                    tid = re.findall(args.prefix + 'T[\d]+', feature_info)[0];
                    transcript_type = re.findall('transcript_biotype "[\w]+"', feature_info)[0].replace("transcript_biotype ", "").replace("\"", "") + " ";

                genes[gid]['transcripts'][tid] = { 'chr' : chrome, 'start' : start, 'end' : end, 'strand' : strand, 'type' : transcript_type, 'cds' : {} };
                num_transcripts += 1;
        core.PWS(core.spacedOut("# Transcripts read:", pad) + str(num_transcripts), outfile);
        core.PWS("# ----------------", outfile);
        # Read transcripts.

    #if not args.no_exons:
        core.PWS("# " + core.getDateTime() + " Reading CDS...", outfile);

        num_cds = 0;
        exon_counts = defaultdict(int);
        for line in reader(args.input, read_mode):
            #print(line);
            line = line_reader(line);
            if "##FASTA" in line[0]:
                break;
            if line[0][0] == "#":
                continue;
            feature_type, chrome, start, end, strand, feature_info = line[2], line[0], line[3], line[4], line[6], line[8];

            if feature_type == "exon":
                if args.source == "maker":
                    #gid = re.findall(args.prefix + '[\d_G]+', feature_info)[0];
                    tid = re.findall(args.prefix + '_T[\d]+_mRNA-[\d]+', feature_info)[0];
                    for gid in genes:
                        if tid in genes[gid]['transcripts']:
                            break; 
                    exon_num = str(len(genes[gid]['transcripts'][tid]['cds']) + 1);
                    eid = tid + "-" + exon_num
                    exon_type = "";
                elif args.source == "ensembl":
                    #print(line);
                    gid = re.findall(args.prefix + 'G[\d]+', feature_info)[0];
                    tid = re.findall(args.prefix + 'T[\d]+', feature_info)[0];
                    eid = re.findall(args.prefix + 'E[\d]+', feature_info);
                    if eid == []:
                        exon_counts[tid] += 1;
                        eid = tid + "-" + str(exon_counts[tid]);
                    else:
                        eid = eid[0];
                    exon_type = re.findall('transcript_biotype "[\w]+"', feature_info)[0].replace("transcript_biotype ", "").replace("\"", "") + " ";
                    
                genes[gid]['transcripts'][tid]['cds'][eid] = { 'chr' : chrome, 'start' : start, 'end' : end, 'strand' : strand, 'type' : exon_type };
                num_cds += 1;
        core.PWS(core.spacedOut("# CDS read:", pad) + str(num_cds), outfile);
        core.PWS("# ----------------", outfile);
        # Read exons.


    core.PWS("# " + core.getDateTime() + " Writing output...", outfile);
    headers = "feature id\tfeature type\tchr\tstart\tend\tstrand\tparent ids\tnum child features";
    outfile.write(headers + "\n");
    # Output headers.

    num_genes, num_transcripts, num_exons = 0,0,0;
    for gid in genes:
        if not args.no_genes:
            outline = [gid, "gene", genes[gid]['chr'], genes[gid]['start'], genes[gid]['end'], genes[gid]['strand'], "NA", str(len(genes[gid]['transcripts']))];
            outfile.write("\t".join(outline) + "\n");
            num_genes += 1;
        # Gene output.

        if not args.no_transcripts or not args.no_exons:
            for tid in genes[gid]['transcripts']:
                cur_transcript = genes[gid]['transcripts'][tid];
                if not args.no_transcripts:
                    outline = [tid, cur_transcript['type'] + "transcript", cur_transcript['chr'], cur_transcript['start'], cur_transcript['end'], cur_transcript['strand'], gid, str(len(cur_transcript['cds']))];
                    outfile.write("\t".join(outline) + "\n");
                    num_transcripts += 1
            # Transcript output.

                if not args.no_exons:
                    cds_starts = { eid : int(cur_transcript['cds'][eid]['start']) for eid in cur_transcript['cds'] };
                    cds_starts = {k: v for k, v in sorted(cds_starts.items(), key=lambda item: item[1])};
                    # Sort the starting coordinates of the CDS.

                    cds_ids = list(cds_starts.keys());
                    if cur_transcript['strand'] == "-":
                        cds_ids.reverse();
                    # Reverse the order if the strand is -.

                    for eid in cds_ids:
                        cur_exon = cur_transcript['cds'][eid];
                        outline = [eid, cur_exon['type'] + "exon", cur_exon['chr'], cur_exon['start'], cur_exon['end'], cur_exon['strand'], gid + ";" + tid, "NA"];
                        outfile.write("\t".join(outline) + "\n");
                        num_exons += 1;
                # Exon output.

    core.PWS(core.spacedOut("# Genes written:", pad) + str(num_genes), outfile);
    core.PWS(core.spacedOut("# Transcripts written:", pad) + str(num_transcripts), outfile);
    core.PWS(core.spacedOut("# Exons written:", pad) + str(num_exons), outfile);
    core.PWS("# Done!", outfile);
    core.PWS("# ----------------", outfile);