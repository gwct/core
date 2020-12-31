############################################################
# Functions for the free-ratio (M1) model from codeml.
# 11.2020
############################################################

import sys, os, re, lib.pamlcore as pcore, lib.pamlseq as pseq, lib.pamltree as ptree

############################################################

def template():
    ctlfile_template = '''seqfile = {infile}
treefile = {treefile}
outfile = {outfile}

model = 2
NSsites = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
aaRatefile = codeml/dat/wag.dat

icode = 0
fix_kappa = 0
kappa = 3
fix_omega = 0
omega = 1

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
getSE = 0
Small_Diff = .5e-6

noisy = 3
verbose = 0
'''
    return(ctlfile_template);

############################################################

def assignTargetClade(tree_str, targs):
# Function that reads a tree into a dictionary and determines the target branch from a set of tip labels

    tinfo, t, root = ptree.treeParse(tree_str);
    # Parse tree into dictionary

    target_found, repl_nodes = False, [];
    for n in tinfo:
        cur_clade = ptree.getClade(n, tinfo);
        # Get the clade that defines each node

        if all(c_nodes in targs for c_nodes in cur_clade):
            repl_nodes.append(n);
            target_found = True;
        # If the clade matches the specified targets, set as a target branch

    if target_found:
        for n in repl_nodes:
            t = t.replace(n, n + " #1");
        t = re.sub('<[\d]+>', '', t) + ";";
    # Add target labels and remove node labels

    return target_found, t;

############################################################

def generate(indir, tree_input, gt_opt, targets, paml_path, outdir, outfile):
    ctlfile_template = template();

    aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
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