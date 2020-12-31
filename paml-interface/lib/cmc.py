############################################################
# Functions for the clade model C (CMC) model from codeml.
# 11.2020
############################################################

import sys, os, re, lib.pamlcore as pcore, lib.pamlseq as pseq, lib.pamltree as ptree

############################################################

def template():
    ctlfile_template = '''seqfile = {infile}
treefile = {treefile}
outfile = {outfile}

model = 3
NSsites = 2

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
ncatG = 3
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

    target_found, target_tree = False, "";
    for n in tinfo:
        cur_clade = ptree.getClade(n, tinfo);
        # Get the clade that defines each node
        if targs == set(cur_clade):
            target_tree = t.replace(n, n + " $1");
            target_found = True;
            break;
        # If the clade matches the specified targets, set as target branch

    if target_found:
        target_tree = re.sub('<[\d]+>', '', target_tree) + ";";
    # Remove node labels

    return target_found, target_tree;

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
        # Read the tree and remove any bootstrap node labels.1
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
        headers = ["file","id","chr","start","end","lnl","k","prop 1", "prop 2", "prop 3", "background clade omega 1", "background clade omega 2", "background clade omega 3", "target clade omega 1", "target clade omega 2", "target clade omega 3"];
    else:
        headers = ["file","lnl","k","prop 1", "prop 2", "prop 3", "background clade omega 1", "background clade omega 2", "background clade omega 3", "target clade omega 1", "target clade omega 2", "target clade omega 3"];
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
                'lnl' : "NA", 'k' : "NA", "prop 1" : "NA", "prop 2" : "NA", "prop 3" : "NA", 
                "background clade omega 1" : "NA", "background clade omega 2" : "NA", "background clade omega 3" : "NA", 
                "target clade omega 1" : "NA", "target clade omega 2" : "NA", "target clade omega 3" : "NA" };
        else:
            gene_info = { 'lnl' : "NA", 'k' : "NA", "prop 1" : "NA", "prop 2" : "NA", "prop 3" : "NA", 
                "background clade omega 1" : "NA", "background clade omega 2" : "NA", "background clade omega 3" : "NA", 
                "target clade omega 1" : "NA", "target clade omega 2" : "NA", "target clade omega 3" : "NA" }
        # Add meta headers to the output dict if metadata is provided

        read_paml_rates = False;
        # Flags to indicate when we've reached certain points in the PAML output

        mono_skipped = False;

        cur_codeml_file = os.path.join(cur_dir, "codeml.out");
        for line in open(cur_codeml_file):
            if line.startswith("lnL"):
                line = list(filter(None, line.strip().split(" ")));
                gene_info['lnl'] = line[4];
                continue;
            # This line indicates that the next line will contain the PAML ancestral relationship labels (ANC..DESC)

            if line.startswith("kappa (ts/tv)"):
                line = line.strip().split(" ");
                #print(line);
                gene_info['k'] = line[-1];
                continue;
            # Get the kappa value for the gene

            if line.startswith("MLEs of dN/dS (w) for site classes (K=3)"):
                read_paml_rates = True;

            if read_paml_rates:
                if line.startswith("proportion"):
                    line = list(filter(None, line.strip().split(" ")));
                    gene_info["prop 1"] = line[1];
                    gene_info["prop 2"] = line[2];
                    gene_info["prop 3"] = line[3];
                    continue;

                if line.startswith("branch type 0"):
                    line = list(filter(None, line.strip().split(" ")));
                    gene_info["background clade omega 1"] = line[3];
                    gene_info["background clade omega 2"] = line[4];
                    gene_info["background clade omega 3"] = line[5];
                    continue;

                if line.startswith("branch type 1"):
                    line = list(filter(None, line.strip().split(" ")));
                    gene_info["target clade omega 1"] = line[3];
                    gene_info["target clade omega 2"] = line[4];
                    gene_info["target clade omega 3"] = line[5];
                    read_paml_rates = False;
                    continue;
                continue;


        gene_outline = [d] + [ gene_info[h] for h in headers if h != "file" ];
        outfile.write(",".join(gene_outline) + "\n");
        # Output all info for current gene

    pcore.PWS("# ----------------", outfile);
    #pcore.PWS(pcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);