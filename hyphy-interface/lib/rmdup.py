############################################################
# Functions for the remove duplicates script in HyPhy.
# 09.2021
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def generate(indir, tree_input, model_file, gt_opt, aln_id_delim, hyphy_path, outdir, logdir, outfile):
    #model_file = os.path.abspath("hyphy-analyses/remove-duplicates/remove-duplicates.bf");

    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    trimmed_treedir = os.path.join(outdir, "trimmed-trees");
    if not os.path.isdir(trimmed_treedir):
        os.system("mkdir " + trimmed_treedir);

    tree_skipped, stop_skipped = 0, 0;
    for aln in aligns:
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
        else:
            tree_file = tree_input;

        if os.path.isfile(tree_file):
            aligns[aln]['tree'] = tree_file;
        # Read the tree and remove any bootstrap node labels.
        # Read the appropriate tree depending on the -tree and -genetree options.

        else:
            outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
            tree_skipped += 1;
            continue;
        # Check the tree file.

        print(aln);
        print(aligns[aln]['aln-file']);
        print(aligns[aln]['tree']);

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

        cur_aln_out = os.path.join(outdir, aln + "-rmdup.fa");
        cur_tree_out = os.path.splitext(os.path.basename(aligns[aln]['tree']));
        cur_tree_out = cur_tree_out[0] + ".trimmed" + cur_tree_out[1];
        cur_tree_out = os.path.join(trimmed_treedir, cur_tree_out);

        #cur_jsonfile = os.path.join(outdir, aln + ".json");
        #cur_outfile = os.path.join(cur_outdir, align + "-out.txt");
        cur_logfile = os.path.join(logdir, aln + ".log");
        # Get the control and output file names

        hyphy_cmd = "hyphy " + model_file + " --msa " + aligns[aln]['aln-file'] + " --tree " +  aligns[aln]['tree'] + " --output " + cur_aln_out + " ENV=\"DATA_FILE_PRINT_FORMAT=9\" &> " + cur_logfile;
        hyphy_cmd += "; tail -n1 " + cur_aln_out + " > " + cur_tree_out;
        hyphy_cmd += "; sed -i '$ d' " + cur_aln_out;
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################
# Keeping parse in case we want to get some stats from the output
'''
def parse(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","dn/ds","lrt","pval"];
    else:
        headers = ["file","branch","dn/ds","lrt","pval"];
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
                "dn/ds" : "NA", "lrt" : "NA", "pval" : "NA" };
        else:
            gene_info = { "dn/ds" : "NA", "lrt" : "NA", "pval" : "NA" };   
        # Initialize the output dictionary for the current branch.

        gene_info["dn/ds"] = str(cur_data["fits"]["Standard MG94"]["Rate Distributions"]["non-synonymous/synonymous rate ratio"]);
        gene_info["lrt"] = str(cur_data["test results"]["non-synonymous/synonymous rate ratio"]["LRT"]);
        gene_info["pval"] = str(cur_data["test results"]["non-synonymous/synonymous rate ratio"]["p-value"]);
        # Retrieve the rate estimates from the json data.

        gene_outline = [f] + [ gene_info[h] for h in headers if h not in ["file"] ];
        outfile.write(",".join(gene_outline) + "\n");
        # Ouput rate estimates for both the gene.

    hpcore.PWS("# ----------------", outfile);
    #hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);
'''