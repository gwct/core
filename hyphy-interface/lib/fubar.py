############################################################
# Functions for the BUSTED model from Hyphy.
# 11.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def generate(indir, tree_input, gt_opt, aln_id_delim, hyphy_path, outdir, logdir, outfile):
    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    for aln in aligns:
        if gt_opt:
            if aln_id_delim:
                tree_dir = os.path.join(tree_input, aln);
                if os.path.isdir(tree_dir):
                    tree_dir_files = os.listdir(tree_dir);
                    tree_file = "".join([ f for f in tree_dir_files if re.findall(aligns[aln]['id'] + '(.*).treefile', f) != [] and "rooted" not in f ]);
                    tree_file = os.path.join(tree_dir, tree_file);
                else:
                    tree_file = False;
            else:
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
        cur_cachefile = os.path.join(outdir, aln + ".cache");
        #cur_outfile = os.path.join(cur_outdir, align + "-out.txt");
        cur_logfile = os.path.join(logdir, aln + ".log");
        # Get the control and output file names

        #hyphy_cmd = "hyphy fubar --alignment " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'] + " --cache " + cur_cachefile + " --output " + cur_jsonfile " &> " + cur_logfile 
        #Call hyphy from conda install
        hyphy_cmd = "/home/ek112884/software/hyphy/hyphy /home/ek112884/software/hyphy/res/TemplateBatchFiles/SelectionAnalyses/FUBAR.bf --alignment " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'] + " --cache " + cur_cachefile + " --output " + cur_jsonfile + " &> " + cur_logfile;
        #Call hyphy from develop branch with FUBAR output file bug fixed
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################

def parse(indir, features, outfile, sitesfile, pad):

    sf = open(sitesfile, "w");

    if features:
        headers = ["file","id","chr","start","end","num ps sites","num ns sites"];
        site_headers = ["file", "id", "site", "site key"];
    else:
        headers = ["file","num ps sites","num ns sites"];
        site_headers = ["file", "site", "site key"];
    outfile.write(",".join(headers) + "\n");
    #sitesfile.write(",".join(site_headers) + "\n");
    sf.write(",".join(site_headers) + "\n");
    # Write the output headers 

    hyphy_files = os.listdir(indir);
    hyphy_files = [ f for f in hyphy_files if f.endswith(".json") ];
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);
    # Read align file names from input directory

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    counter = 0;
    for f in hyphy_files:
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
                            "num ps sites" : 0, "ps sites" : [], "num ns sites" : 0, "ns sites" : [] };
        else:
            gene_info = { "num ps sites" : 0, "ps sites" : [], "num ns sites" : 0, "ns sites" : [] };   
        # Initialize the output dictionary for the current branch.

        #gene_info["dn/ds"] = str(cur_data["fits"]["Standard MG94"]["Rate Distributions"]["non-synonymous/synonymous rate ratio"]);

        site_mat = cur_data["MLE"]["content"]["0"];
        site_ind = 1;
        for site in site_mat:
            if site[3] > 0.9:
                gene_info["num ns sites"] += 1;
                gene_info["ns sites"].append(site_ind);
            # Negative selection

            if site[4] > 0.9:
                gene_info["num ps sites"] += 1;
                gene_info["ps sites"].append(site_ind);

                site_str = str(site_ind);
                if features:
                    site_info = { 'file' : f, 'id' : fid, 'site' : site_str, 'site key' : fid + ":" + site_str };
                else:
                    site_info = { 'file' : f, 'site' : site_str, 'site key' : f + ":" + site_str };
                site_outline = [ site_info[h] for h in site_headers ];
                #sitesfile.write(",".join(site_outline) + "\n");
                sf.write(",".join(site_outline) + "\n");
            # Positive selection

            site_ind += 1;
        # Retrieve the positively and negatively selected sites from the json data.

        gene_outline = [f] + [ str(gene_info[h]) for h in headers if h not in ["file"] ];
        #gene_outline = [ ";".join(c) if type(c) == list else c for c in gene_outline ];
        outfile.write(",".join(gene_outline) + "\n");
        # Ouput rate estimates for both the gene.

    hpcore.PWS("# ----------------", outfile);
    #hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);
    sf.close()
