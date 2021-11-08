############################################################
# General R functions for figure development,, 03.20
#############################################################

bartheme <- function () {  
# My standard theme for most plots.
  theme_classic() %+replace% 
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=16), 
          axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black",angle=90), 
          axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0),color="black"),
          axis.line=element_line(colour='#595959',size=0.75),
          axis.ticks=element_line(colour="#595959",size = 1),
          axis.ticks.length=unit(0.2,"cm"),
          legend.position="right",
          legend.key.width = unit(0.75,  unit = "cm"),
          legend.spacing.x = unit(0.25, 'cm'),
          legend.title = element_blank(),
          legend.text=element_text(size=12),
          plot.title = element_text(hjust=0.5, size=20),
          plot.margin = unit(c(1,1,1,1), "cm")
    )
}

#############################################################

dottheme <- function () {
# A theme for dotplots.
  theme_linedraw() %+replace% 
    theme(axis.text=element_text(size=8), 
          axis.title=element_text(size=16), 
          axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0),color="black",angle=90), 
          axis.title.x=element_text(margin=margin(t=0,r=0,b=10,l=0),color="black"),
          #axis.line=element_blank(),#element_line(colour='#595959',size=0.75),
          axis.ticks=element_line(colour="#595959", size=0.25),
          axis.ticks.length=unit(0.1,"cm"),
          panel.grid.major = element_line(color="#999999"),
          panel.grid.minor = element_line(color="#bfbfbf"),
          legend.position="right",
          legend.key.width = unit(0.75,  unit = "cm"),
          legend.spacing.x = unit(0.25, 'cm'),
          legend.title = element_text(size=12),
          legend.text=element_text(size=10),
          plot.title = element_text(hjust=0.5, size=20),
          plot.margin = unit(c(1,1,1,1), "cm")
    )
}

############################################################

corecol <- function (pal="default", numcol=4, offset=0, info=FALSE) {
# Custom color palettes.
    palette_list = c("default", "trek", "trekdark", "wilke")

    if(!pal %in% palette_list){
        cat("# * CORECOL: Requested palette not in palette list.\n")
        cat("# * CORECOL: Requested:    ", pal, "\n")
        cat("# * CORECOL: Palette list: ", palette_list, "\n")
        cat("# * CORECOL: Setting palette to 'default'\n")
        pal = "default"
    }

    if(pal=="default"){
        col = c("#db6d00", "#004949", "#006ddb", "#920000", "#490092", "#6cb6ff", "#24ff24", "#fdb4da", "#ffff6d", "#009292", "#924900", "#000000")
    }else if(pal=="trek"){
        col = c("#4A508A", "#D89000", "#C44040", "#4d3d7e", "#6c465b", "#78736f", "#000000")
    }
    else if(pal=="trekdark"){
        col = c("#62121f", "#7d5811", "#174159", "#4d3d7e", "#6c465b", "#78736f", "#000000")
    }else if(pal=="wilke"){
        col = c("#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7", "#000000")
    }

    if(numcol > length(col)){
        cat("# * CORECOL: Not enough colors in selected palette:\n")
        cat("# * CORECOL: numcol     = ", numcol, "\n")
        cat("# * CORECOL: pal        = ", pal, "\n")
        cat("# * CORECOL: pal length = ", length(col), "\n")
    }

    return_col = col

    offset_counter = offset
    while(offset_counter>0){
        return_col = c(return_col, return_col[1])
        # Copies first element at end

        return_col = return_col[2:(length(return_col))]
        # Removes first element

        offset_counter = offset_counter - 1
    }

    return_col = head(return_col, numcol)

    if(info){
        cat("\n# CORECOL INFO\n")
        cat("# Palette list:               ", palette_list, "\n")
        cat("# Requested palette:          ", pal, "\n")
        cat("# Original colors:            ", col, "\n")
        cat("# Num colors in palette:      ", length(col), "\n")
        cat("# Num colors requested:       ", numcol, "\n")
        cat("# Offset:                     ", offset, "\n")
        cat("# Returned colors:            ", return_col, "\n")
        cat("# ------------------\n")
    }

    return(return_col)
}

############################################################

nodecheck <- function (tree, tree_type="object", xmax=1) {
# A function to check ggtree's internal node labels.

  if(!tree_type %in% c("object", "string", "file")){
      cat("\n# nodecheck: Invalid tree_type specified:    ", tree_type, "\n")
      cat("\n# nodecheck: tree_type must be one of:       object, string, file\n")
      cat("# ------------------\n")     
  }

  if(tree_type == "string"){
      tree = read.tree(text=tree_str)
  }else if(tree_type == "file"){
      tree = read.tree(tree_file)
  }

  node_test = ggtree(tree, size=2, ladderize=F) +
    ggplot2::xlim(0, xmax) +
    geom_tiplab(color="#333333", fontface='italic', size=5) +
    geom_text(aes(label=node), hjust=-.3, vjust=-.3, color="#ff6db6") +
    geom_nodepoint(color="#666666", alpha=0.85, size=4) +
    ggtitle("NODECHECK (if no tree is visible, try changing xmax from default of 1)")
  print(node_test)
}