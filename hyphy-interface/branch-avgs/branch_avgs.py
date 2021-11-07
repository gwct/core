#!/usr/bin/env python3
#############################################################################
# Computes branch averages for rates from a directory of
# hyphy results as csv files.
# Modified from Gregg Thomas's script (31 August 2021) to parse ES, EN, S, and N instead of dN/dS (EK)
# 10.29.2021: Reformatted into a more project like script (GT)
#############################################################################

import sys
import os
import multiprocessing as mp
import lib.core as CORE
import lib.params as params
import lib.opt_parse as OP
import lib.branches as BRANCHES

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# branch_avgs version " + globs['version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.

    print("#");
    print("# " + "=" * 125);
    #print(CORE.welcome());
    #if "-h" not in sys.argv:
    #    print("            Degeneracy annotation of transcripts\n");
    # A welcome banner.

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ##########################

    step = "Reading files in input directory";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    globs['csv-files'] = [ f for f in os.listdir(globs['csv-rate-dir']) if f.endswith(".csv") ];
    # Read the filtered genes to skip them for averaging

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(globs['csv-files'])) + " files read");
    # Status update

    ##########################

    if globs['filter-file']:
        step = "Reading filter file";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Status update

        globs['filter-files'] = [ line.replace("\"", "").replace(".json", ".csv") for line in open(globs['filter-file']) ];
        # Read the filtered genes to skip them for averaging

        globs['csv-files'] = list( set(globs['csv-files']) - set(globs['filter-files']) );
        # Remove the files to filter from the list of input files

        step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(globs['filter-files'])) + " loci will be excluded");
        # Status update
    
    ##########################

    if globs['subset-file']:
        step = "Reading subset file";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Status update

        globs['subset-files'] = [ line.replace("\n","-mafft-cds.filter.csv") for line in open(globs['subset-file']) ];
        # Read in a subset of genes to include in analysis

        globs['csv-files'] = list( set(globs['csv-files']) & set(globs['subset-files']) );
        # Only take files in both lists

        step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(globs['subset-files'])) + " loci will be included");
        print(len(globs['csv-files']));
        # Status update

    ########################## 

    step = "Reading species tree information";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    branches = {};
    # The main dictionary for the tree info

    first = True;
    branch_num = 1;
    for line in open(globs['tree-info-file']):
        line = line.strip().replace("\"", "").split(",");
        # Parse the line of the file

        if first:
            globs['headers'] = line;        
            first = False;
            continue;
        # If it's the first line, get the headers and skip

        clade = line[globs['clade-index']]
        globs['branches'][clade] = {};
        branches[clade] = {};
        # Get the clade and add it as a key to the main dictionary

        for i in range(len(globs['headers'])):
        # Add each header value into the dictionary for the current clade
            if i == globs['clade-index']:
                continue;
            # The clade will actually be the key, so we don't need to add it into the dictionary

            globs['branches'][clade][globs['headers'][i]] = line[i];
            branches[clade][globs['headers'][i]] = line[i];
            # Add the value

        globs['branches'][clade]['branch.num'] = str(branch_num);
        branches[clade]['branch.num'] = str(branch_num);
        branch_num += 1;

    globs['num-branches'] = len(globs['branches']);
    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(globs['num-branches']) + " branches read");
    # Status update

    # Read in the tree info
    ##########################

    step = "Adding new headers";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update
    new_headers = ["ES.sum", "EN.sum", "S.sum", "N.sum", "num.genes.full", "num.genes.partial", "num.genes.no.clade"];
    # The purpose of this script is to add some new columns to the tre csv, as defined here

    for clade in globs['branches']:
    # For every clade in branches, add the new headers into the branches dictionary

        for nh in new_headers:
            if globs['branches'][clade]["node.type"] == "ROOT":
                globs['branches'][clade][nh] = "NA";
                branches[clade][nh] = "NA";
            else:
                globs['branches'][clade][nh] = 0;
                branches[clade][nh] = 0;
        # Add in the new headers and initialize at 0, except for the root node

    globs['headers'] += new_headers;
    # Combine the original and new headers for the output

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(new_headers)) + " headers added");
    # Status update

    # Compile the headers for the new columns
    ##########################

    step = "Summing values over branches";
    step_start_time = CORE.report_step(globs, step, False, "Processed 0 / " + str(globs['num-branches']) + " branches...", full_update=True);
    # Status update

    # branch_num = 1;
    # for branch in branches:
    #     if branch != "Kadarsanomys_sodyi_MZB-Sample":
    #         continue;

    #     result = BRANCHES.branchSum((globs['branches'][branch], branch, globs['csv-rate-dir'], globs['filter-files'], globs['subset-files']));
    #     new_branches[result[0]] = result[1];
    #     #branches.update(result);
    #     cur_branch_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(branch_num) + " / " + str(globs['num-branches']) + " branches...", full_update=True);
    #     branch_num += 1;
    ## Serial version for debugging

    with mp.Pool(processes=globs['num-procs']) as pool:
        branch_num = 1;
        new_branches = {};
        for result in pool.imap_unordered(BRANCHES.branchSum, ((globs['branches'][branch], branch, globs['csv-rate-dir'], globs['csv-files'], globs['filter-files'], globs['subset-files']) for branch in branches)):
            #print(result[0], result[1]);
            
            new_branches[result[0]] = result[1];
            #branches.update(result);
            cur_branch_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(branch_num) + " / " + str(globs['num-branches']) + " branches...", full_update=True);
            branch_num += 1;
    ## Parallel version

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=True);
    # Status update

    # Get number of subs for each branch from each gene
    ##########################

    step = "Writing out new tree table";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    avg_headers = ["avg.ES", "avg.EN", "avg.S", "avg.N", "dS", "dN", "dNdS"];
    globs['headers'] += avg_headers;

    for clade in new_branches:
    # For every clade in branches, add the new headers into the branches dictionary

        for nh in avg_headers:
            if new_branches[clade]["node.type"] == "ROOT":
                new_branches[clade][nh] = "NA";
            else:
                new_branches[clade][nh] = 0;
        # Add in the new headers and initialize at 0, except for the root node

    with open(globs['output-file'], "w") as outfile:
        outfile.write(",".join(globs['headers']) + "\n")
        # Open the output file and write the headers, which now contain the new columns

        for branch in new_branches:
            if new_branches[branch]["node.type"] != "ROOT" and (new_branches[branch]['num.genes.full'] + new_branches[branch]['num.genes.partial']) != 0:
                new_branches[branch]['avg.ES'] = new_branches[branch]['ES.sum']  / (new_branches[branch]['num.genes.full'] + new_branches[branch]['num.genes.partial']);
                new_branches[branch]['avg.EN'] = new_branches[branch]['EN.sum']  / (new_branches[branch]['num.genes.full'] + new_branches[branch]['num.genes.partial']);
                new_branches[branch]['avg.S'] = new_branches[branch]['S.sum']  / (new_branches[branch]['num.genes.full'] + new_branches[branch]['num.genes.partial']);
                new_branches[branch]['avg.N'] = new_branches[branch]['N.sum']  / (new_branches[branch]['num.genes.full'] + new_branches[branch]['num.genes.partial']);
            # If the branch is not the root and appears in some genes, then compute the averages.

                new_branches[branch]['dS'] = new_branches[branch]['S.sum'] / new_branches[branch]['ES.sum']
                new_branches[branch]['dN'] = new_branches[branch]['N.sum'] / new_branches[branch]['EN.sum']
                new_branches[branch]['dNdS'] = new_branches[branch]['dN'] / new_branches[branch]['dS']

                # new_branches[branch]['avg.dS'] = new_branches[branch]['avg.S'] / new_branches[branch]['avg.ES']
                # new_branches[branch]['avg.dN'] = new_branches[branch]['avg.N'] / new_branches[branch]['avg.EN']
                # new_branches[branch]['avg.dNdS'] = new_branches[branch]['avg.dN'] / new_branches[branch]['avg.dS']
            # After averaging ES, EN, S, and N, use them to calculate dS, dN, and dN/dS across all alignments

            outline = [];
            for header in globs['headers']:
                if header == 'clade':
                    outline.append(branch);
                else:
                    outline.append(str(new_branches[branch][header]));
            
            outline = ",".join(outline);
            outfile.write(outline + "\n");
            # Compile the output line and write it out to the file.

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    # Output the tree table with the new columns
    ##########################

    CORE.endProg(globs);

#############################################################################

