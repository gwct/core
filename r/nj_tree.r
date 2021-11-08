########################################################################################
#A little script to make a quick NJ tree in R.
#
#Dependencies: ape
#
#Gregg Thomas, Summer 2015
########################################################################################

args=(commandArgs(TRUE))
library(ape)

sdm_distances = read.table(file=args[1])
sdm_matrix = as.matrix(sdm_distances)
sdm_tree = nj(sdm_matrix)
write.tree(sdm_tree, args[2])
