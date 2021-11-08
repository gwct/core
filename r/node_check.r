############################################################
# Checks ggtree node labels for an input tree
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
library(ggtree)

args=(commandArgs(TRUE))

#tree_file = "C:/bin/hamster-genome/data/mol-evol/mus-4spec-iqtree-concat-rooted.tre"
tree_file = args[1]
tree = read.tree(tree_file)

node_test = ggtree(tree, size=2, ladderize=F) +
  ggplot2::xlim(0, 0.05) +
  geom_tiplab(color="#333333", fontface='italic', size=5) +
  geom_text(aes(label=node), hjust=-.3, vjust=-.3, color="#ff6db6") +
  geom_nodepoint(color="#666666", alpha=0.85, size=4)
print(node_test) 