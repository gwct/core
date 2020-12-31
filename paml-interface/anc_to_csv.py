#!/usr/bin/python
############################################################
# Reads results from codeml ancestral reconstructions.(rst file).
############################################################

import sys, os, argparse, re, lib.pamlcore as pcore, lib.pamltree as ptree

############################################################
# Options

parser = argparse.ArgumentParser(description="Parse PAML codeml output for a tip and target branches");
parser.add_argument("-i", dest="input", help="Directory containing subdirectories of codeml runs.", default=False);
#parser.add_argument("-m", dest="model", help="The PAML model that was used to generate the files in -i. Default: m1", default="m1");
parser.add_argument("-d", dest="meta", help="A file containing metadata. Tab delimited with columns: id, feature type, chromosome, start coord, end coord, strand. Directories in -i must be formatted <id>-*", default=False);
parser.add_argument("-o", dest="output", help="An output .csv file.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit(" * Error 1: Please provide a valid input directory (-i).");

# if args.model not in ["m1", "m2", "cmc"]:
#     sys.exit(" * Error 2: Model (-m) must be one of: m1, m2, cmc");

if args.meta and not os.path.isfile(args.meta):
    sys.exit(" * Error 3: Cannot find meta data file: " + args.meta);

if not args.output:
    sys.exit(" * Error 2: Please provide the name of an output file (-o).")

if os.path.isfile(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output file (-o) already exists! Explicity specify --overwrite to overwrite it.");

pad = 25;

with open(args.output, "w") as outfile:
    pcore.runTime("# codeml m1 output parser", outfile);
    pcore.PWS("# IO OPTIONS", outfile);
    pcore.PWS(pcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    #pcore.PWS(pcore.spacedOut("# PAML model:", pad) + args.model, outfile);
    if args.meta:
        pcore.PWS(pcore.spacedOut("# Metadata file:", pad) + args.meta, outfile);
    pcore.PWS(pcore.spacedOut("# Output file:", pad) + args.output, outfile);
    if args.overwrite:
        pcore.PWS(pcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous output file.", outfile);
    pcore.PWS("# ----------------", outfile);

    pcore.PWS("# " + pcore.getDateTime() + " Reading BLOSUM62 matrix to calculate avg. scores per branch", outfile);
    blosum62 = {};
    first = True;
    for line in open("lib/BLOSUM62"):
        if line[0] == "#":
            continue;
        if first:
            headers = list(filter(None, line.strip().split(" ")));
            first = False;
            continue;
        
        line = list(filter(None, line.strip().split(" ")));
        aa1 = line[0];
        scores = line[1:];
        for i in range(len(scores)):
            aa2 = headers[i];
            blosum62[aa1 + "-" + aa2] = float(scores[i]);
    pcore.PWS("# ----------------", outfile);
    # Read the BLOSUM62 matrix from a text file into a dictionary.

    features = False;
    if args.meta:
        pcore.PWS("# " + pcore.getDateTime() + " Reading metadata file: " + args.meta, outfile);
        features = pcore.readMeta(args.meta);
        pcore.PWS(pcore.spacedOut("# Features read:        ", pad) + str(len(features)), outfile);               
        pcore.PWS("# ----------------", outfile);
    # Read the feature metadata.

    if features:
        headers = ["file","id","chr","start","end","clade", "num subs", "avg sub score"];
    else:
        headers = ["file","clade","num subs", "avg sub score"];
    outfile.write(",".join(headers) + "\n");
    # Write the output headers 

    num_unfinished = 0;
    # A count of the number of unfinished codeml runs as determined by incomplete output files

    num_nonmono = 0;
    # A count of the number of trees where the target clade is not monophyletic.

    codeml_dirs = os.listdir(args.input);
    num_dirs = len(codeml_dirs);
    num_dirs_str = str(num_dirs);
    num_dirs_len = len(num_dirs_str);
    # Read align file names from input directory

    counter = 0;
    for d in os.listdir(args.input):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_dirs_len:
                counter_str = "0" + counter_str;
            print ("> " + pcore.getDateTime() + " " + counter_str + " / " + num_dirs_str);
        counter += 1;
        # Track progress

        #print(d);
        cur_dir = os.path.join(args.input, d);
        # Get the current gene ID and codeml directory.

        if features:
            if "-" in d:
                fid = d.split("-")[0];
            elif "." in d:
                fid = d.split(".")[0];
            else:
                fid = d;
            cur_feature = features[fid];
            # Look up transcript/gene info for current gene to save in output

            gene_outline = [d, fid, cur_feature['chrome'], cur_feature['start'], cur_feature['end']];
        else:
            gene_outline = [d];
        # Get the gene output line if metadata features are provided.

        read_paml_ancs = True;
        read_paml_tree = False;
        read_branches = False;
        # Flags to indicate when we've reached certain points in the PAML output

        paml_ancs = {};
        # This will be the main tree info dictionary, akin to the treeParse dictionary
        # <PAML node label> : [<branch length>, <ancestral node>, <node type>, <original tip node label>]

        cur_codeml_file = os.path.join(cur_dir, "rst");
        #print(cur_codeml_file);
        for line in open(cur_codeml_file):
            if read_paml_ancs and ".." in line:
            # The first line with .. in it indicates ancestral relationships of nodes in the tree in PAML's format.

                nodes = list(filter(None, line.strip().split(" ")));
                # Get the relationships into a list

                for node in nodes:
                    anc, n = node.split("..");
                    paml_ancs[n] = ['', anc, '', ''];
                # Read through every PAML ANC..DESC relationship and add it to the main dictionary.

                for node in nodes:
                    anc, n = node.split("..");
                    if anc not in paml_ancs:
                        paml_ancs[anc] = ['NA', "NA", 'root', ''];
                        root = anc;
                # The trifurcating 'root' node will not be included in the loop above (since it has no ancestor). Add it to the
                # dictionary as 'root' here.

                read_paml_ancs = False;
                # Need to indicate that we don't want to run this block on the subsequent lines

            if line.startswith("tree with node labels for Rod Page's TreeView"):
                read_paml_tree = True;
                continue;
            # This line indicates that the next line will contain the tree with PAML tip labels

            if read_paml_tree:
            # Read the PAML tree to get the orignal tip labels associated with PAML's labels
                for node in paml_ancs:
                # Parse every node in the main dictionary

                    if node + "_" in line:
                    # If the node is found in the tree string with an underscore it is a tip node.
                        paml_ancs[node][2] = 'tip'
                        # Add this node as a tip in the main dictionary

                        orig_label = re.findall(node + "_" + "[\w,\)]+", line)[0][:-1];
                        paml_ancs[node][3] = orig_label[orig_label.index("_")+1:];
                        # Add the original tip label into the paml_ancs dictionary.

                    elif paml_ancs[node][2] != 'root':
                        paml_ancs[node][2] = 'internal'
                    # If the node isn't a tip or the trifurcating 'root', add it as internal with no original label saved.

                max_clade_node, max_clade_count = "", 0;
                for n in paml_ancs:
                # For every internal node we now assign clades.

                    if paml_ancs[n][2] != 'tip':
                        cur_clade_paml = ptree.getClade(n, paml_ancs);
                        cur_clade_orig = sorted([ paml_ancs[tip][3] for tip in cur_clade_paml ]);
                        paml_ancs[n][3] = cur_clade_orig;

                        if paml_ancs[n][2] != 'root' and len(cur_clade_paml) > max_clade_count:
                            max_clade_node, max_clade_count = n, len(cur_clade_paml);
                        # Get the longest clade other than the "root" clade

                for n in paml_ancs:
                    if paml_ancs[n][2] == 'tip':
                        paml_ancs[n][3] = [paml_ancs[n][3]];
                # Convert the tip labels into lists    

                root_clade = sorted(list(set(paml_ancs[root][3]) - set(paml_ancs[max_clade_node][3])));
                paml_ancs[root][3] = root_clade;
                # Subtract out the longest clade from the 'root' clade to get the remaining clade.

                read_paml_tree = False;
                continue;

            if line.startswith("List of extant and reconstructed sequences"):
            # This line indicates we're done reading branches. This has to come before the actual branch reading block so we don't skip it.
                for outdict in [branch_1_info, branch_2_info]:
                    avg_sub_score = 0.0;
                    if outdict['num subs'] != 0:
                        avg_sub_score = outdict['sum sub score'] / outdict['num subs'];
                    # Calculate the average BLOSUM score for the branch.

                    branch_outline = gene_outline + [ outdict['clade'], str(outdict['num subs']), str(avg_sub_score) ];
                    outfile.write(",".join(branch_outline) + "\n");
                # The last branch read is output here.
                
                read_branches = False;
                break;
                # This should be the last line we need to read, so exit the loop here.

            if line.startswith("Branch"):
            # For every line starting with "Branch" we have to 1) Output the previous branch and 2) Initialize the new branch.
                
                if read_branches:
                    for outdict in [branch_1_info, branch_2_info]:
                        avg_sub_score = 0.0;
                        if outdict['num subs'] != 0:
                            avg_sub_score = outdict['sum sub score'] / outdict['num subs'];
                        # Calculate the average BLOSUM score for the branch.

                        branch_outline = gene_outline + [ outdict['clade'], str(outdict['num subs']), str(avg_sub_score) ];
                        outfile.write(",".join(branch_outline) + "\n");
                # This block outputs the previous branch if read branches is true (which won't be the case for the first time through when we haven't
                # read any branches yet).

                read_branches = True;
                # Set read branches to true so we know to read the next lines.

                line = list(filter(None, line.strip().split(" ")));
                anc, node = line[2].split("..");
                assert paml_ancs[node][1] == anc, "\nMismatching ancestral node!\n" + d + "\n" + paml_ancs + "\n" + anc + " " + node;
                # Get the nodes on joining the current branch and confirm their relationship

                clade = paml_ancs[node][3];
                anti_clade = sorted([ paml_ancs[n][3][0] for n in paml_ancs if paml_ancs[n][2] == 'tip' and paml_ancs[n][3][0] not in clade ]);
                # Get the clades on both sides of the branch.

                branch_1_info = { 'clade' : ";".join(clade), "num subs" : 0.0, "sum sub score" : 0.0 };
                branch_2_info = { 'clade' : ";".join(anti_clade), "num subs" : 0.0, "sum sub score" : 0.0 };
                # Initialize the two clades output dicts.
             
                continue;

            if read_branches and line != "\n":
            # Every non-blank line that defines a substitution.

                line = list(filter(None, line.strip().split(" ")));
                anc_allele = line[2].replace("(","").replace(")","");
                der_allele = line[6].replace("(","").replace(")","");
                # Parse the alleles int he substitution.

                if anc_allele == der_allele or "-" in [anc_allele, der_allele]:
                    continue;
                # If there is no non-synonymous substitution or its an indel, skip.

                branch_1_info['num subs'] += 1.0;
                branch_1_info['sum sub score'] += blosum62[anc_allele + "-" + der_allele];
                # Get the score for one direction.

                branch_2_info['num subs'] += 1.0;
                branch_2_info['sum sub score'] += blosum62[der_allele + "-" + anc_allele];
                # Get the score for the other direction.