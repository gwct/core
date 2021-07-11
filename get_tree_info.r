############################################################
# Gets the information from a tree to start a data file, 06.20
#############################################################

library(optparse)
library(phytools)
library(ape)

#############################################################

treeToDF <- function (tree, tree_type="object") {
# A function to extract info from a tree and save it as a data frame.
# Info includes: node labels, clades, node types, and branch lengths.

  if(!tree_type %in% c("object", "string", "file")){
    cat("\n# treeToDF: Invalid tree_type specified:    ", tree_type, "\n")
    cat("\n# treeToDF: tree_type must be one of:       object, string, file\n")
    cat("# ------------------\n")     
  }
  # Make sure a valid tree type has been passed.
  
  if(tree_type == "string"){
    tree = read.tree(text=tree_str)
  }else if(tree_type == "file"){
    tree = read.tree(tree_file)
  }
  # If the tree type is a string or a file, read it as a phy object here
  
  clades = data.frame("node"=c(), "clade"=c(), "node.type"=c())
  # Data frame of clades for each node in the tree. A clade is the tip labels that 
  # descend from the node.
  
  done = c()
  for(n in tree[["edge"]][,1]){
  # Go through every internal node in the tree (including the root)
    
    if(!n %in% done){
      done = c(done, n)
      # If we haven't seen this node before, add it to the list of ones we've seen
      
      node_type = "internal"
      # Assign default node type as internal
      
      desc = getDescendants(tree, n)
      # Get all descendants of this node
      
      clade = c()
      for(d in desc){
        if(d <= length(tree[["tip.label"]])){
          clade = c(clade, tree[["tip.label"]][d])
        }
      }
      # The tips are assigned node labels first so the node number matches the index in the
      # tip.label vector. This converts between node number and tip label.
      
      if(identical(sort(clade), sort(tree[["tip.label"]]))){
        node_type = "ROOT"
      }
      # If the length of the clade equals the number of tips, assign this as the root
      
      clades = rbind(clades, data.frame("node"=n, "clade"=paste(sort(clade), collapse=";"), "node.type"=node_type))
      # Add the current clade to the clades data frame
    }
  }
  
  for(i in 1:length(tree[["tip.label"]])){
    clades = rbind(clades, data.frame("node"=i, "clade"=tree[["tip.label"]][i], "node.type"="tip"))
  }
  # Add the tip labels as their own clades
  
  
  blens = data.frame("node"=c(), "branch.length"=c())
  for(i in 1:length(tree[["edge.length"]])){
    blens = rbind(blens, data.frame("node"=tree[["edge"]][i,2], "branch.length"=tree[["edge.length"]][i]))
  }
  # Extract the branch lengths for all nodes except the root
  
  tree_info = merge(clades, blens, by="node", all=T)
  # Combine the clade and branch length data frames
  
  tree_info$label = NA
  for(i in 1:nrow(tree_info)){
    #row = tree_info[i,]
    #print(row)
    
    if(tree_info[i,]$node.type == "tip"){
      tree_info[i,]$label = tree_info[i,]$clade
    }else{
      tree_info[i,]$label = tree$node.label[tree_info[i,]$node - Ntip(tree)]
    }
  }

  return(tree_info)
}

option_list = list(
  make_option(c("-t", "--treefile"), action="store", default=FALSE, type='character',
              help="A file containing a newick tree. Required."),
  make_option(c("-o", "--outfile"), action="store", default=FALSE, type='character',
              help="A .csv file to write the data frame."),
  make_option(c("-w", "--overwrite"), action="store_true", default=FALSE,
              help="If the output file already exists, this must be set to overwrite it."),
  make_option(c("-q", "--quiet"), action="store_true", default=FALSE,
              help="Set to not print the tree info to the screen.")
)
opt = parse_args(OptionParser(option_list=option_list))

if(opt[["treefile"]] == "FALSE" && opt[["outfile"]] == "FALSE" && opt[["help"]] == FALSE){
  warning(" * WARNING 1: Provide arguments. Use -h for options.")
}
if(opt[["treefile"]] == "FALSE" || !file.exists(opt[["treefile"]])){
  warning(" * WARNING 2: A file with a Newick formatted tree must be provided with -t.")
}
if(file.exists(opt[["outfile"]]) && !opt[["overwrite"]]){
  stop(" * ERROR 1: Output file (-o) already exists! Explicity specify --overwrite to overwrite it.")
}

read_tinfo = FALSE
if(opt[["treefile"]] != "FALSE"){
  cat(as.character(Sys.time()), " | Reading tree in file: ", opt[["treefile"]], "\n")
  tree = read.tree(opt[["treefile"]])
  tinfo = treeToDF(tree)
  read_tinfo = TRUE
}

if(read_tinfo && !opt[["quiet"]]){
  cat("-------------------------------------\n")
  print(tinfo)
  cat("-------------------------------------\n")
}

if(opt[["outfile"]] != "FALSE"){
  cat(as.character(Sys.time()), " | Writing tree info to file: ", opt[["outfile"]], "\n")
  write.csv(tinfo, file=opt[["outfile"]], row.names=F)
}








