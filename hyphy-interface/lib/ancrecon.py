############################################################
# Functions for the free-ratio (MG94) model from Hyphy.
# 11.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

def generate(indir, model_file, hyphy_path, outdir, logdir, outfile):
    for f in os.listdir(indir):
        if not f.endswith(".json"):
            continue;

        infilename = os.path.join(indir, f);
        outfilename = os.path.join(outdir, f.replace(".json", ".ancestors.json"));
        logfilename = os.path.join(logdir, f.replace(".json", ".log"));

        hyphy_cmd = "hyphy " + model_file + " --fit " + infilename + " --output " + outfilename + " &> " + logfilename;
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

############################################################

def parse(indir, features, outfile, pad):

    out_mode = "clades";
    print("out mode: " + out_mode);
    # clades or splits

    outdir = os.path.join(indir, "csv");
    if not os.path.isdir(outdir):
        os.system("mkdir " + outdir);

    if features:
        if out_mode == "splits":
            headers = ["file","id","chr","start","end","branch","split1","split2","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
        elif out_mode == "clades":
            headers = ["file","id","chr","start","end","branch","clade","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
    else:
        if out_mode == "splits":
            headers = ["file","branch","split1","split2","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
        elif out_mode == "clades":
            headers = ["file","branch","clade","dn","ds","dn/ds","nonsynonymous bl","synonymous bl"];
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

        gene_outfile = os.path.join(outdir, f.replace(".json", ".csv"));
        with open(gene_outfile, "w") as goutfile:
            goutfile.write(",".join(headers) + "\n");

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

            splits = defaultdict(list);
            for n in tinfo:
                if tinfo[n][2] == 'root':
                    continue;
                if tinfo[n][2] == 'tip':
                    node_label = n;
                else:
                    node_label = tinfo[n][3];

                split1 = tp.getClade(n, tinfo);
                split1.sort();
                #print(cur_clade);
                split2 = [ n2 for n2 in tinfo if tinfo[n2][2] == 'tip' and n2 not in split1 ];
                split2.sort();

                splits[node_label].append(";".join(split1));
                splits[node_label].append(";".join(split2));
                # For unrooted trees there is no directionality, so get clades from both sides of each branch.
            # For each node/branch in the tree, save the tip species that make up the clade.

            for node in splits:
                if features:
                    if out_mode == "splits":
                        branch_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                            "branch" : node, "split1" : "NA", "split2" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };
                    elif out_mode == "clades":
                        branch_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                            "branch" : node, "clade" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };                        
                else:
                    if out_mode == "splits":
                        branch_info = { "branch" : node ,"split1" : "NA", "split2" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };   
                    elif out_mode == "clades":
                        branch_info = { "branch" : node ,"clades" : "NA", "dn" : "NA" ,"ds" : "NA" ,"dn/ds" : "NA", "nonsynonymous bl" : "NA", "synonymous bl" : "NA" };   
                # Initialize the output dictionary for the current branch.

                branch_info["dn/ds"] = str(cur_data["branch attributes"]["0"][node]["Confidence Intervals"]["MLE"]);
                branch_info["dn"] = str(cur_data["branch attributes"]["0"][node]["dN"]);
                branch_info["ds"] = str(cur_data["branch attributes"]["0"][node]["dS"]);
                branch_info["nonsynonymous bl"] = str(cur_data["branch attributes"]["0"][node]["nonsynonymous"]);
                branch_info["synonymous bl"] = str(cur_data["branch attributes"]["0"][node]["synonymous"]);
                # Retrieve the rate estimates from the json data.

                if out_mode == "splits":
                    branch_info["split1"] = splits[node][0];
                    branch_info["split2"] = splits[node][1];

                    branch_outline = [f] + [ branch_info[h] for h in headers if h not in ["file"] ];
                    outfile.write(",".join(branch_outline) + "\n");
                    goutfile.write(",".join(branch_outline) + "\n");

                if out_mode == "clades":
                    for clade in splits[node]:
                        branch_info["clade"] = clade;
                        branch_outline = [f] + [ branch_info[h] for h in headers if h not in ["file"] ];
                        outfile.write(",".join(branch_outline) + "\n");
                        goutfile.write(",".join(branch_outline) + "\n");
                # Ouput rate estimates for both clades of the current branch.
    hpcore.PWS("# ----------------", outfile);
    hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);