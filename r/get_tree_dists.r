############################################################
# Gets tree distance between all trees in two files, 10.21
#############################################################

library(optparse)
library(phytools)
library(ape)
library(phangorn)
library(stringr)

#############################################################

getDists <- function(t1, t2, dist_type){
  
  cat(as.character(Sys.time()), " | Reading tree files\n")
  trees1 = readLines(t1, warn=F)
  trees2 = readLines(t2, warn=F)
  # Read the tree files
  
  dist_col = paste(dist_type, ".dist", sep="")
  # A string with the current distance type
  
  out_data = data.frame("tree1"=c(), "tree2"=c(), "dist_col"=c())
  # Initialize the data frame with a generic "dist_col" that will get replaced later
  
  cat(as.character(Sys.time()), " | Calculating distances\n")
  for(tree_str1 in trees1){
    tree1 = read.tree(text=tree_str1)
    # Read the first tree as a phylo object
    
    for(tree_str2 in trees2){
      tree2 = read.tree(text=tree_str2)
      # Read the second tree as a phylo object
      
      if(dist_type == "rf"){
        cur_dist = RF.dist(tree1, tree2)
      }else if(dist_type == "wrf"){
        cur_dist = wRF.dist(tree1, tree2)
      }else if(dist_type == "spr"){
        cur_dist = SPR.dist(tree1, tree2)
      }else if(dist_type == "kf"){
        cur_dist = KF.dist(tree1, tree2)
      }else if(dist_type == "path"){
        cur_dist = path.dist(tree1, tree2)
      }
      # Calculate the distance between the two trees based on the input option
      
      out_data = rbind(out_data, data.frame("tree1"=tree_str1, "tree2"=tree_str2, "dist_col"=cur_dist))
      # Add the current distance to the data frame
    }
    ## End second tree loop
  }
  ## End first tree loop
  
  names(out_data) = c("tree1", "tree2", dist_col)
  # Replace the generic distance column header with the specific one
  
  return(out_data)
}

#############################################################

option_list = list(
  make_option(c("-t", "--treefile1"), action="store", default=FALSE, type='character',
              help="A file containing one or more newick trees. Required."),
  make_option(c("-r", "--treefile2"), action="store", default=FALSE, type='character',
              help="A file containing one or more newick trees. Required."),
  make_option(c("-d", "--disttype"), action="store", default=FALSE, type='character',
              help="One of 'rf', 'wrf', 'spr', 'kf', or 'path'. Required."),
  make_option(c("-o", "--outfile"), action="store", default=FALSE, type='character',
              help="A file to write out the RFs"),
  make_option(c("-w", "--overwrite"), action="store_true", default=FALSE,
              help="If the output file already exists, this must be set to overwrite it."),
  make_option(c("-q", "--quiet"), action="store_true", default=FALSE,
              help="Set to not print the tree info to the screen.")
)
opt = parse_args(OptionParser(option_list=option_list))
# Input options

opt[["treefile1"]] = "C:/Users/grt814/Desktop/phodopus_4spec_iqtree.cf.tree"
opt[["treefile2"]] = "C:/Users/grt814/Desktop/loci.treefile"
opt[["disttype"]] = "rf"
#opt[["outfile"]] = "test.csv"


if(opt[["treefile1"]] == "FALSE" && opt[["treefile2"]] == "FALSE" && opt[["disttype"]] == "FALSE" && opt[["outfile"]] == "FALSE" && opt[["help"]] == FALSE){
  warning(" * WARNING 1: Provide arguments. Use -h for options.")
}
if(opt[["treefile1"]] == "FALSE" || !file.exists(opt[["treefile1"]])){
  warning(" * WARNING 2: A file with a Newick formatted tree must be provided with -t1.")
}
if(opt[["treefile2"]] == "FALSE" || !file.exists(opt[["treefile2"]])){
  warning(" * WARNING 3: A file with a Newick formatted tree must be provided with -t2.")
}
if(!opt[["disttype"]] %in% c("rf", "wrf", "spr", "kf", "path")){
  warning(" * WARNING 4: -d must be one of 'rf', 'wrf', 'spr', 'kf', or 'path'.")
}
if(file.exists(opt[["outfile"]]) && !opt[["overwrite"]]){
  stop(" * ERROR 1: Output file (-o) already exists! Explicity specify --overwrite to overwrite it.")
}
# Input option checking

dists = getDists(opt[["treefile1"]], opt[["treefile2"]], opt[["disttype"]])
# Call the distance function

dist_col = paste(opt[["disttype"]], ".dist", sep="")
# The column header for the current distance option

hist(dists[[dist_col]])
# A histogram of distances

if(opt[["outfile"]] != "FALSE"){
  cat(as.character(Sys.time()), " | Writing distances: ", opt[["outfile"]], "\n")
  write.csv(dists, file=opt[["outfile"]], row.names=F)
}
# Output to a file, if specified

#############################################################
