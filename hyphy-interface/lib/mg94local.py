############################################################
# Functions for the free-ratio (MG94) model from Hyphy.
# 11.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def generate(indir, tree_input, gt_opt, hyphy_path, outdir, logdir, outfile):
    model_file = os.path.abspath("hyphy-analyses/FitMG94/FitMG94.bf");

    aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    for aln in aligns:
        if gt_opt:
            tree_file = os.path.join(tree_input, aln, aln + ".treefile");
        else:
            tree_file = tree_input;

        if os.path.isfile(tree_file):
            aligns[aln]['tree'] = tree_file;
        # Read the tree and remove any bootstrap node labels.1
    # Read the appropriate tree depending on the -tree and -genetree options.

    tree_skipped, stop_skipped = 0, 0;
    for aln in aligns:
        if not aligns[aln]['tree']:
            outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
            tree_skipped += 1;
            continue;
        # Check the tree file.          

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

        hyphy_cmd = "hyphy " + model_file + " --alignment " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'] + " --type local --output " + cur_jsonfile + " &> " + cur_logfile 
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################

def parse(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","branch","clade","dn","ds","dn/ds"];
    else:
        headers = ["file","branch","clade","dn","ds","dn/ds"];
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

        cur_tree = cur_data['input']['trees']['0'];
        tinfo, tree, root = tp.treeParse(cur_tree);
        # Read the tree from the json data.

        clades = defaultdict(list);
        for n in tinfo:
            if tinfo[n][2] == 'root':
                continue;
            if tinfo[n][2] == 'tip':
                node_label = n;
            else:
                node_label = tinfo[n][3];

            cur_clade = tp.getClade(n, tinfo);
            cur_anti_clade = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in cur_clade ];
            clades[node_label].append(";".join(cur_clade));
            clades[node_label].append(";".join(cur_anti_clade));
            # For unrooted trees there is no directionality, so get clades from both sides of each branch.
        # For each node/branch in the tree, save the tip species that make up the clade.

        for node in clades:
            if features:
                branch_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                    "branch" : node, "clade" : "NA" ,"dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA" };
            else:
                branch_info = { "branch" : node ,"clade" : "NA" ,"dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA" };   
            # Initialize the output dictionary for the current branch.

            branch_info["dn/ds"] = str(cur_data["branch attributes"]["0"][node]["Confidence Intervals"]["MLE"]);
            branch_info["dn"] = str(cur_data["branch attributes"]["0"][node]["nonsynonymous"]);
            branch_info["ds"] = str(cur_data["branch attributes"]["0"][node]["synonymous"]);
            # Retrieve the rate estimates from the json data.

            for clade in clades[node]:
                branch_info["clade"] = clade;
                branch_outline = [f] + [ branch_info[h] for h in headers if h not in ["file"] ];
                outfile.write(",".join(branch_outline) + "\n");
            # Ouput rate estimates for both clades of the current branch.
    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);