############################################################
# Functions for the aBSREL model from Hyphy.
# 12.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def assignTargetClade(tree_file, targets):
# Function that reads a tree into a dictionary and determines the target branch from a set of tip labels

    tree_str = open(tree_file, "r").read().strip();
    tinfo, t, root = tp.treeParse(tree_str);
    # Parse tree into dictionary
    
    #print(tinfo);
    #print(t);
    target_label = "NA";
    target_found = True;

    max_split1 = [];
    max_split2 = [];
    max_node = "";
    # Variables to store the information for the maximal subset clade

    for node in tinfo:
        split1 = tp.getClade(node, tinfo);
        split1.sort();
        split2 = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in split1 ];
        split2.sort();

        splits = [split1, split2];
        # Get the splits from the current branch in the gene's gene tree

        for s in range(len(splits)):
        # For each split, we want to check whether it is a subset of the current branch

            if set(splits[s]) == set(targets) or set(splits[s]) <= set(targets):
            # Check if the current split is an equivalent set or subset of the targets

                if len(splits[s]) > len(max_split1):
                # Check if the current split that is a subset of the targets contains more species
                # than the current maximal subset. If so, it becomes the maximal subset.

                    max_split1 = splits[s];
                    max_node = node;
                    if s == 0:
                        max_split2 = splits[1];
                    else:
                        max_split2 = splits[0];
                    # Assign the maximal subset to max_split1 and the other split to max_split2
    # Check every node in the tree for target species

    if max_split1 == [] or any(spec in max_split2 for spec in targets):
        target_found = False;
    # If no species from the current branch are found then this clade doesn't exist in this gene OR
    # If species from the branch are present in the second split (not the maximal subset), then this is a case of
    # discordance and the clade truly doesn't exist in this gene as a monophyly. Do not increment sums.

    else:
        target_label = max_node;
        print(target_label)
    # Otherwise save the node as the target label

        # if tinfo[max_node][2] == 'tip':
        #     target_label = max_node;
        # else:
        #     target_label = "<target #1>"
        #     t = t.replace(max_node, target_label);

    return t, target_label, target_found;

############################################################

def generate(indir, tree_input, gt_opt, aln_id_delim, target_clade, hyphy_path, outdir, target_tree_dir, logdir, outfile):
    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False, "target" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False, "target" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    #tree_skipped, stop_skipped = 0, 0;
    #for aln in aligns:
    #    if not aligns[aln]['tree']:
    #        outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
    #        tree_skipped += 1;
    #        continue;
        # Check the tree file.
   
    #tree_skipped, stop_skipped = 0, 0;
    #if not os.path.exists(tree_input):    
    #    for aln in aligns:    
    #        outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
    #        tree_skipped += 1;
    #        continue;
        # Check the tree file.

    if target_clade and not gt_opt:
        target_tree, target_label, target_found,  = assignTargetClade(tree_input, target_clade);
        # Assign the target label for the input tree
        if not target_found:
            sys.exit(" * ERROR: The target clade does not exist in the provided tree.");
        # If the target branch doesn't exist in the single input tree, exit
        else:
            tree_input = os.path.basename(tree_input);
            tree_input = os.path.join(target_tree_dir, tree_input);
        # Write the labeled input tree for the run
    # This block handles the target branch identification for a single input tree

    no_target_trees = 0;
    # A count of the number of trees that don't have the target branch, for gene tree input

    #tree_skipped, stop_skipped = 0, 0;
    #for aln in aligns:
    #    if not aligns[aln]['tree']:
    #        outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
    #        tree_skipped += 1;
    #        continue;
    #    # Check the tree file.

    tree_skipped, stop_skipped = 0, 0;
    for aln in aligns:
        print(aln)
        if gt_opt:
            if aln_id_delim:
                tree_dir = os.path.join(tree_input, aln);
                if os.path.isdir(tree_dir):
                    tree_dir_files = os.listdir(tree_dir);
                    tree_file = "".join([ f for f in tree_dir_files if re.findall(aligns[aln]['id'] + '(.*).treefile', f) != [] and "rooted" not in f ]);
                    if not tree_file:
                        continue;
                    tree_file = os.path.join(tree_dir, tree_file);
                else:
                    tree_file = False;
            # If we need to split the input alignment directory to get the tree file name
            else:
                tree_file = os.path.join(tree_input, aln, aln + ".treefile");
            # Get the tree file name
            
            if os.path.isfile(tree_file):
                aligns[aln]['tree'] = tree_file;
            # Assign the current tree file to the alignment
            if not aligns[aln]['tree']:
                print(aligns[aln]['tree'])
                outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
                tree_skipped += 1;
                continue;
            # Check the tree file.

            if target_clade:
                target_tree, target_label, target_found = assignTargetClade(tree_file, target_clade);
                # Identify the target node of the current tree

                if not target_found:
                    no_target_trees += 1;
                # Count the tree if the target clade isn't found
                else:
                    tree_file = os.path.basename(tree_file).replace(".treefile", "-labeled.treefile");
                    tree_file = os.path.join(target_tree_dir, tree_file);
                    aligns[aln]['target'] = target_label;
                    with open(tree_file, "w") as tout:
                        tout.write(target_tree);
                # Re-write the labeled tree file if the target clade is found
            else:
                target_label = False;
        else:
            tree_file = tree_input;
        # Gene tree target lineage labeling block

        if os.path.isfile(tree_file):
            aligns[aln]['tree'] = tree_file;
            aligns[aln]['target'] = target_label;
            # Assign the current tree file and target label to the alignment
	    # Read the appropriate tree depending on the -tree and -genetree options.

	    #tree_skipped, stop_skipped = 0, 0;
	    #for aln in aligns:
	    #    if not aligns[aln]['tree']:
	    #        outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
	    #        tree_skipped += 1;
	    #        continue;
	    #    # Check the tree file.          

            seq_dict = hpseq.fastaGetDict(aligns[aln]['aln-file']);
            prem_stop_flag = False
            for title in seq_dict:
                prem_stop, new_seq = hpseq.premStopCheck(seq_dict[title], allowlastcodon=True, rmlast=True);
                if prem_stop:
                    prem_stop_flag = True;
                    seq_dict[title] = new_seq;
                # Read the sequence and check for premature stop codons.

                if prem_stop_flag:
                    outfile.write(" # Premature stop found. Skipping: " + aln + "\n");
                    stop_skipped += 1;
                    continue;
                # Check the alignment for premature stop codons (which PAML hangs on)

        # cur_outdir = os.path.join(outdir, aln);
        # if not os.path.isdir(cur_outdir):
        #     os.system("mkdir " + cur_outdir);
        # Make the output directory for this alignment

        cur_jsonfile = os.path.join(outdir, aln + ".json");
        #cur_outfile = os.path.join(cur_outdir, align + "-out.txt");
        cur_logfile = os.path.join(logdir, aln + ".log");
        # Get the control and output file names

        hyphy_cmd = "hyphy absrel --alignment " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'];
        if target_clade:
            if not aligns[aln]['target']:
                continue;
            hyphy_cmd += " --branches " + aligns[aln]['target'];
        # Add the target node if specified and found
        hyphy_cmd += " --output " + cur_jsonfile + " &> " + cur_logfile 
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    if target_clade:
        hpcore.PWS("# Num skipped because target clade not found     : " + str(no_target_trees), outfile);
    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);
    # Report some stats

############################################################

def parse(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","num ps branches", "ps pvals"];
    else:
        headers = ["file","branches","num ps branches","ps pvals"];
    outfile.write(",".join(headers) + "\n");
    # Write the output headers 

    hyphy_files = os.listdir(indir);
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);
    # Read align file names from input directory

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    counter = 0;
    for f in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
        counter += 1;
        # Track progress

        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".json"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current output file.

        #print(f);
        if features:
            if "-" in f:
                fid = f.split("-")[0];
            elif "." in f:
                fid = f.split(".")[0];
            else:
                fid = f;
            cur_feature = features[fid];
            # Look up transcript/gene info for current gene to save in output, if metadata is provided.

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the Hyphy json file.

        #print(cur_data)

        if features:
            gene_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                "num ps branches" : 0, "ps pvals" : [] };
        else:
            gene_info = { "branches" : [], "num ps branches" : 0, "ps pvals" : [] };   
        # Initialize the output dictionary for the current branch.

        #gene_info["dn/ds"] = str(cur_data["fits"]["Standard MG94"]["Rate Distributions"]["non-synonymous/synonymous rate ratio"]);
        for node in cur_data["branch attributes"]["0"]:
            if float(cur_data["branch attributes"]["0"][node]["Corrected P-value"]) < 0.01:
                gene_info["branches"].append(str(node)); 
                gene_info["num ps branches"] += 1;
                gene_info["ps pvals"].append(str(cur_data["branch attributes"]["0"][node]["Corrected P-value"]));
        # Retrieve the rate estimates from the json data.

        gene_info["branches"] = ";".join(gene_info["branches"]);
        gene_info["ps pvals"] = ";".join(gene_info["ps pvals"]);
        gene_info["num ps branches"] = str(gene_info["num ps branches"]);
        gene_outline = [f] + [ gene_info[h] for h in headers if h not in ["file"] ];
        outfile.write(",".join(gene_outline) + "\n");
        # Ouput rate estimates for both the gene.

    hpcore.PWS("# ----------------", outfile);
    #hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);






























































'''
def generate(indir, tree_input, gt_opt, targets, aln_id_delim, hyphy_path, outdir, logdir, outfile):

    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    for aln in aligns:
        if gt_opt:
            tree_file = os.path.join(tree_input, aln, aln + "-rooted.treefile");
        else:
            tree_file = tree_input;

        if os.path.isfile(tree_file):
            cur_tree = open(tree_file, "r").read().strip();
            cur_tree = re.sub("\)[\de.-]+:", "):", cur_tree);
            aligns[aln]['tree'] = cur_tree;
        # Read the tree and remove any bootstrap node labels.
    # Read the appropriate tree depending on the -tree and -genetree options.

    tree_skipped, target_skipped, stop_skipped = 0, 0, 0;
    for aln in aligns:
        if not aligns[aln]['tree']:
            outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
            tree_skipped += 1;
            continue;
        # Check the tree file.          

        target_found, cur_tree = assignTargetClade(aligns[aln]['tree'], targets);
        if not target_found:
            outfile.write(" # Target clade not found. Skipping: " + aligns[aln]['aln-file'] + "\n");
            target_skipped += 1;
            continue;
        # Assign the target lineages to the current tree.

        seq_dict = pseq.fastaGetDict(aligns[aln]['aln-file']);
        prem_stop_flag = False
        for title in seq_dict:
            prem_stop, new_seq = pseq.premStopCheck(seq_dict[title], allowlastcodon=True, rmlast=True);
            if prem_stop:
                prem_stop_flag = True;
            seq_dict[title] = new_seq;
        # Read the sequence and check for premature stop codons.

        if prem_stop_flag:
            outfile.write(" # Premature stop found. Skipping: " + aln + "\n");
            stop_skipped += 1;
            continue;
        # Check the alignment for premature stop codons (which PAML hangs on)

        cur_outdir = os.path.join(outdir, aln);
        if not os.path.isdir(cur_outdir):
            os.system("mkdir " + cur_outdir);
        # Make the output directory for this alignment

        new_treefile = os.path.join(cur_outdir, "codeml.tre");
        with open(new_treefile, "w") as treefile:
            treefile.write(cur_tree);
        # Make the tree file for this alignment

        new_seqfile = os.path.join(cur_outdir, "codeml.fa");
        with open(new_seqfile, "w") as seqfile:
            for title in seq_dict:
                seqfile.write(title + "\n");
                seqfile.write(seq_dict[title] + "\n");
        # Write the sequences for this alignment

        cur_ctlfile = os.path.join(cur_outdir, "codeml.ctl");
        cur_outfile = os.path.join(cur_outdir, "codeml.out");
        #cur_logfile = os.path.join(cur_outdir, base_input + "-codeml.log");
        # Get the control and output file names

        with open(cur_ctlfile, "w") as ctlfile:
            ctlfile.write(ctlfile_template.format(infile=new_seqfile, treefile=new_treefile, outfile=cur_outfile));
        # Write the control file

        codeml_cmd = "cd " + cur_outdir + "; " + paml_path + " codeml.ctl";
        outfile.write(codeml_cmd + "\n");
        # Construct and write the codeml command

    pcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    pcore.PWS("# Num skipped because target clade not found  : " + str(target_skipped), outfile);
    pcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################

def parse(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","lnl","k","dn sum","ds sum","dn avg","ds avg","clade","clade dn/ds","clade dn","clade ds"];
    else:
        headers = ["file","lnl","k","dn sum","ds sum","dn avg","ds avg","clade","clade dn/ds","clade dn","clade ds"];
    outfile.write(",".join(headers) + "\n");
    # Write the output headers 

    num_unfinished = 0;
    # A count of the number of unfinished codeml runs as determined by incomplete output files

    num_nonmono = 0;
    # A count of the number of trees where the target clade is not monophyletic.

    codeml_dirs = os.listdir(indir);
    num_dirs = len(codeml_dirs);
    num_dirs_str = str(num_dirs);
    num_dirs_len = len(num_dirs_str);
    # Read align file names from input directory

    counter = 0;
    for d in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_dirs_len:
                counter_str = "0" + counter_str;
            print ("> " + pcore.getDateTime() + " " + counter_str + " / " + num_dirs_str);
        counter += 1;
        # Track progress

        cur_dir = os.path.join(indir, d);
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

            gene_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                'lnl' : "NA", 'k' : "NA", 'dnds' : "NA", 'dn sum' : "NA", 'ds sum' : "NA", 'dn avg' : "NA", 'ds avg' : "NA" };
        else:
            gene_info = { 'lnl' : "NA", 'k' : "NA", 'dnds' : "NA", 'dn sum' : "NA", 'ds sum' : "NA", 'dn avg' : "NA", 'ds avg' : "NA" }
        # Add meta headers to the output dict if metadata is provided

        branch_info = {};
        # Initialize the output dictionary for the current gene

        read_paml_ancs = False;
        read_paml_tree = False;
        read_orig_tree = False;
        # Flags to indicate when we've reached certain points in the PAML output

        num_branches = 0.0;
        # Count the number of branches in the tree to calculate per branch averages

        paml_ancs = {};
        # This will be the main tree info dictionary, akin to the treeParse dictionary
        # <PAML node label> : [<branch length>, <ancestral node>, <node type>, <original tip node label>]

        mono_skipped = False;

        cur_codeml_file = os.path.join(cur_dir, "codeml.out");
        for line in open(cur_codeml_file):
            if line.startswith("lnL"):
                line = list(filter(None, line.strip().split(" ")));
                gene_info['lnl'] = line[4];

                read_paml_ancs = True;
                continue;
            # This line indicates that the next line will contain the PAML ancestral relationship labels (ANC..DESC)

            if read_paml_ancs:
            # PAML ancestral relationship labels line.

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


            if line.startswith("tree length ="):
                read_paml_tree = True;
                continue;
            # This line indicates that the next line will contain the tree with PAML tip labels

            if read_paml_tree and line != "\n":
            # Line with tree with PAML tip labels

                for node in paml_ancs:
                # Parse every node in the main dictionary

                    if node + ": " in line:
                    # If the node has a branch length in the tree string, do the following

                        paml_ancs[node][2] = 'tip'
                        # Add this node as a tip in the main dictionary

                        #node_ind = line.index(node + ": ");

                        bl_match = re.findall(node + ': [\d.]+\)', line);
                        if bl_match == []:
                            bl_match = re.findall(node + ': [\d.]+,', line);
                        # Retrieve the branch length for this node based on pattern "<NODE> :<BL>)" OR "<NODE> :<BL>,"
                        
                        cur_bl = bl_match[0];
                        cur_bl = cur_bl[cur_bl.index(" ")+1:-1];
                        paml_ancs[node][0] = cur_bl;
                        # Parse the branch length retrieved and add it to the dictionary

                    elif paml_ancs[node][2] != 'root':
                        paml_ancs[node][2] = 'internal'
                    # If the node isn't a tip or the trifurcating 'root', add it as internal with no branch length saved

                paml_label_order = [ l[:l.index(":")] for l in re.findall('[\d]+: ', line) ];
                # Get all the tip node labels in the order they appear in the tree string

                read_paml_tree = False;
                read_orig_tree = True;
                continue;
                # Indicate that we don't want to run this block again, and that the next line contains the tree with 
                # the original tip labels
                
            if read_orig_tree and line != "\n":
            # Line with the tree with the original tip labels
                #print(line);
                orig_label_order = [ l[:l.index(":")] for l in re.findall('[\w]+: ', line) ];
                # Get the tip node labels in the order they appear in the tree string

                assert len(paml_label_order) == len(orig_label_order), "\nUnequal number of labels: " + d + "\nPAML: " + ",".join(paml_label_order) + "\nOrig:" + ",".join(orig_label_order)
                # Make sure we have the same number of labels in the PAML and original trees

                for n in range(len(paml_label_order)):
                    node = paml_label_order[n];
                    paml_ancs[node][3] = orig_label_order[n];
                # Add the original label into the label field of the main dictionary

                max_clade_node, max_clade_count = "", 0;
                for n in paml_ancs:
                    if paml_ancs[n][2] != 'tip':
                        cur_clade_paml = ptree.getClade(n, paml_ancs);
                        cur_clade_orig = sorted([ paml_ancs[tip][3] for tip in cur_clade_paml ]);
                        paml_ancs[n][3] = cur_clade_orig;

                        if paml_ancs[n][2] != 'root' and len(cur_clade_paml) > max_clade_count:
                            max_clade_node, max_clade_count = n, len(cur_clade_paml);
                        # Get the longest clade other than the "root" clade
                # Get the clades for each node

                for n in paml_ancs:
                    if paml_ancs[n][2] == 'tip':
                        paml_ancs[n][3] = [paml_ancs[n][3]];
                # Convert the tip labels into lists    

                # clade_len_counts = [ len(paml_ancs[n][3]) for n in paml_ancs ];
                # if clade_len_counts.count(max_clade_count) != 1:
                #     print("More than one max clade length found: ");
                #     print(tid);
                #     print(paml_ancs);
                #     print(max_clade_node);
                #     sys.exit(1);
                #assert clade_len_counts.count(max_clade_count) == 1, "\nMore than one max clade length found: " + gid + "\n" + paml_ancs;

                root_clade = sorted(list(set(paml_ancs[root][3]) - set(paml_ancs[max_clade_node][3])));
                paml_ancs[root][3] = root_clade;
                # Subtract out the longest clade from the 'root' clade to get the remaining clade.

                read_orig_tree = False;
                # Indicate that we don't need to run this block again
                continue;

            if line.startswith("kappa (ts/tv)"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['k'] = line[-1];
                continue;
            # Get the kappa value for the gene

            if line.count("..") == 1 and "check convergence" not in line:
            # These lines report the brance specific results

                line = list(filter(None, line.strip().split(" ")));
                # Convert the line into a list

                nodes = line[0].split("..")

                node = nodes[1];
                anc = nodes[0];

                dnds, dn, ds = line[4], line[5], line[6];
                # Parse relevant info from line

                cur_clade = ";".join(paml_ancs[node][3]);
                anti_clade = ";".join(sorted([ s for s in orig_label_order if s not in paml_ancs[node][3] ]));
                branch_info[cur_clade] = { 'dnds' : dnds, 'dn' : dn, 'ds' : ds };
                branch_info[anti_clade] = { 'dnds' : dnds, 'dn' : dn, 'ds' : ds };
                # Since unrooted trees lack directionality, get both clades on both sides of the branch

                # if node == max_clade_node and anc == root:
                #     anc = max_clade_node;
                #     node = root;
                #     cur_clade = ";".join(paml_ancs[node][3]);
                #     branch_info[cur_clade] = { 'dnds' : dnds, 'dn' : dn, 'ds' : ds };
                # For the branch between the 'root' and its opposite, max clade node, reverse the directonality and output that clade as well

                num_branches += 1;
                continue;

            if line.startswith("tree length for dN:"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['dn sum'] = line[-1];
                continue;
            # Sum of dN for all branches

            if line.startswith("tree length for dS:"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['ds sum'] = line[-1];
                continue;
            # Sum of dS for all branches

        try:
            gene_info['dn avg'] = str(float(gene_info['dn sum']) / num_branches);
            gene_info['ds avg'] = str(float(gene_info['ds sum']) / num_branches);   
        except:
            #print(gid);
            if not mono_skipped:
                num_unfinished += 1;
        # Try to calculate average dN and dS for all branches. If this fails, codeml likely didn't finish

        gene_outline = [d] + [ gene_info[h] for h in headers if h not in ["file","clade","clade dn/ds","clade dn","clade ds"] ];

        # gene_outline = [ gene_info['gid'], gene_info['mtid'], gene_info['mchr'], gene_info['mstart'], gene_info['mend'],
        #             gene_info['ptid'], gene_info['pchr'], gene_info['pscaff'], gene_info['pstart'], gene_info['pend'],
        #             gene_info['lnl'], gene_info['k'], 
        #             gene_info['dn-sum'], gene_info['ds-sum'], gene_info['dn-avg'], gene_info['ds-avg'] ];

        for branch in branch_info:
            branch_outline = [ branch, branch_info[branch]['dnds'], branch_info[branch]['dn'], branch_info[branch]['ds'] ];
            outline = gene_outline + branch_outline;
            outfile.write(",".join(outline) + "\n");
        # Output all info for each branch for current gene

    pcore.PWS("# ----------------", outfile);
    pcore.PWS(pcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);
'''
