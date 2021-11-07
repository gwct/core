#############################################################################
# Function for determining whether a branch exists in a gene tree
#############################################################################

import os

#############################################################################

def branchSum(branch_sum_item):
# Retrieve rates/counts for each branch in the species tree for each gene, if that branch exists
# in the gene

    cur_branch_dict, cur_branch, rates_dir, csv_files, filter_files, subset_files = branch_sum_item;

    # print(cur_branch);
    #core.PWS("# " + core.getDateTime() + " Starting branch " + cur_branch);

    if cur_branch_dict["node.type"] == "ROOT":
            return [cur_branch, cur_branch_dict];
    # Skip the root node because it doesn't have an associated branch

    branch_list = cur_branch.split(";");
    # Make a list of the species descendant from the current branch

    num_genes = len(csv_files);
    num_genes_str = str(num_genes);
    # The number of branches and genes, for progress updates, str here to minimize calls to str

    gene_counter = 1;
    # The gene counter for progress updates, the branch counter is also converted to string here to 
    # minimize calls to str

    for f in csv_files:
    # For the current branch, go through every gene and get the maximal subset clade

        if f in filter_files:
            continue;
        # Skip the genes filtered based on dS

        if subset_files and f not in subset_files:
            continue;
        # Skip the genes not in the subset of genes we're interested in

        # if gene_counter == 1 or gene_counter % 100 == 0:
        #     print("# Branch : " + str(cur_branch_num) + " / " + num_branches + " -> gene: " + str(gene_counter) + " / " + num_genes_str);
        gene_counter += 1;
        # Progress updates

        # if gene_counter > 10:
        #     continue;

        infile = os.path.join(rates_dir, f);
        # Compile path to current file

        max_split1 = [];
        max_split2 = [];
        max_line = [];
        # Variables to store the information for the maximal subset clade

        first = True;
        for line in open(infile):
            if first:
                first = False;
                continue;
            # Skip the header line in the file

            line = line.strip().split(",");
            split1 = line[2].split(";");
            split2 = line[3].split(";");
            splits = [split1, split2];
            # Get the splits from the current branch in the gene's gene tree

            for s in range(len(splits)):
            # For each split, we want to check whether it is a subset of the current branch

                if set(splits[s]) == set(branch_list) or set(splits[s]) <= set(branch_list):
                # Check if the current split is an equivalent set or subset of the current branch

                    if len(splits[s]) > len(max_split1):
                    # Check if the current split that is a subset of the branch contains more species
                    # than the current maximal subset. If so, it becomes the maximal subset.

                        max_split1 = splits[s];
                        max_line = line;
                        if s == 0:
                            max_split2 = splits[1];
                        else:
                            max_split2 = splits[0];
                        # Assign the maximal subset to max_split1 and the other split to max_split2

        if max_split1 == [] or any(spec in max_split2 for spec in branch_list):
            cur_branch_dict['num.genes.no.clade'] += 1;
        # If no species from the current branch are found then this clade doesn't exist in this gene OR
        # If species from the branch are present in the second split (not the maximal subset), then this is a case of
        # discordance and the clade truly doesn't exist in this gene as a monophyly. Do not increment sums.

        elif len(branch_list) == len(max_split1):
            cur_ES = float(max_line[4]);
            cur_EN = float(max_line[5]);
            cur_S = float(max_line[6]);
            cur_N = float(max_line[7]);
            # Parse the rates from the line corresponding to the maximal split

            #if cur_ds_bl == 1e-10:
            #    cur_branch_dict['num.genes.no.ds'] += 1;
            #    cur_ds, cur_ds_bl, cur_dnds = 0, 0, 0;
            #elif cur_dn_bl == 1e-10:
            #    cur_branch_dict['num.genes.no.dn'] += 1;
            #    cur_dn, cur_dn_bl, cur_dnds = 0, 0, 0;

            cur_branch_dict['num.genes.full'] += 1;
            cur_branch_dict['ES.sum'] += cur_ES;
            cur_branch_dict['EN.sum'] += cur_EN;
            cur_branch_dict['S.sum'] += cur_S;
            cur_branch_dict['N.sum'] += cur_N;
        # If the maximal subset of the current branch is equal to the current branch, increment sums for each new column.

        else:
            cur_ES = float(max_line[4]);
            cur_EN = float(max_line[5]);
            cur_S = float(max_line[6]);
            cur_N = float(max_line[7]);
            # Parse the rates from the line corresponding to the maximal split

            if len(branch_list) == 1:
                print(cur_branch);
                print(max_split1);
                print(len(max_split1));
                print(max_split2);
                print(len(max_split2));
                sys.exit();
            # This shouldn't happen ... it would mean a tip was found in the opposite split of the maximal clade (which should just be the tip)

            #if cur_ds_bl == 1e-10:
            #    cur_branch_dict['num.genes.no.ds'] += 1;
            #    cur_ds, cur_ds_bl, cur_dnds = 0, 0, 0;
            #elif cur_dn_bl == 1e-10:
            #    cur_branch_dict['num.genes.no.dn'] += 1;
            #    cur_dn, cur_dn_bl, cur_dnds = 0, 0, 0;
            # If there are no substitutions of a given type, the branch length will be 0, which I think Hyphy represents as 1e-10
            # Count these here and set appropriate rates to 0

            cur_branch_dict['num.genes.partial'] += 1;
            cur_branch_dict['ES.sum'] += cur_ES;
            cur_branch_dict['EN.sum'] += cur_EN;
            cur_branch_dict['S.sum'] += cur_S;
            cur_branch_dict['N.sum'] += cur_N;
        # If there are no species from the current branch in the second split (not the maximal subset), then this 
        # is a case of missing data and we can still count the branch as existing in this gene tree. Increment
        # the sums for each new column.

    #core.PWS("# " + core.getDateTime() + " Finishing branch " + cur_branch);
    #print(cur_branch_dict)

    return [cur_branch, cur_branch_dict];

#############################################################################