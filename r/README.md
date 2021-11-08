# core/r/

R scripts for tasks related to handling phylogenetic trees as well as some figure design functions.

Almost all of these scripts are written in R 4.0+ (https://www.r-project.org/).

| File | Description | 
| ------ | ----------- |
| design.r | Scripts that generate commands to run on many files as well as SLURM submission scripts. Includes generators for a coding sequence alignment pipeline. |
| get_tree_dists.r | Takes two tree files as input and computes tree distance (e.g. Robinson-Foulds distance) between all pairs of trees. |
| get_tree_info.r | Reads a Newick formatted tree and converts it into a .csv table with one node per row. Nodes are labeled as ape and ggtree expect them, so this is useful to link the tree to outside data. |
| nj_tree.r | Calculates a neighbor-joining tree from an input distance matrix. |
| node_check.r | Displays a tree with nodes labeled as ape and ggtree expect them. |
