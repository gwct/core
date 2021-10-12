#PURPOSE: Perform multiple test correction on aBSREL output p-values
#NOTE: aBSREL does some kind of multiple test correction already, but this is just among branches within a locus; because we are doing different aBSREL runs for each locus we also need to do a correction across all loci
#NOTE: In theory,this should work for any branch-site tests in HyPhy, but I have only tested with the aBSREL output (EK 2021-10-11)

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: name of HyPhy run (ex: rtm-gt3-absrel)")
}
this_run<-args[1]

#Read in data and get vector of all P-values
mydata<-read.csv(paste0(this_run,".csv"), header=TRUE, comment.char="#")
#print(head(mydata))
myps<-as.character(mydata$ps.pvals)
print(paste("Total number of loci:",length(myps)))
print(paste("Number of loci with at least one branch under positive selection:", length(which(myps!=""))))
#print(myps[1:10])
split_ps<-unlist(sapply(myps, function(x) unlist(strsplit(x, ";"))))
print(paste("Number of significant uncorrected P-values across all branches and loci:", length(split_ps)))
print(split_ps[1:10])

#Get a vector of names for p-value vector; this will keep track of which p-value belongs to which locus and branch
mynames<-c()
for(i in 1:nrow(mydata)){
	counter<-0
	mybranches<-as.character(mydata$branches)[i]
	split_branches<-unlist(strsplit(mybranches, ";"))
	#print(mybranches)
	#print(split_branches)
	#print(length(split_branches))
	#print(counter)
	while(counter < length(split_branches)){
		locus<-gsub("-.*", "", mydata$file[i])
		branch<-split_branches[counter+1]
		mynames<-c(mynames, paste(locus,branch,sep="."))
		counter<-counter+1
	}
	
}
print(length(mynames))
names(split_ps)<-mynames
print(split_ps[1:20])
#print(tail(split_ps))

#Perform FDR correction
fdr.ps<-p.adjust(split_ps, method="fdr")
fdr.ps.sig<-fdr.ps[which(fdr.ps < 0.05)]
print(paste("Length of FDR-corrected P-value vector:", length(fdr.ps)))
print(paste("Number of significant P-values (P<0.05) after FDR correction:", length(fdr.ps.sig)))
print(head(fdr.ps))
fdr.ps_df<-as.data.frame(cbind(names(fdr.ps.sig), fdr.ps.sig))
colnames(fdr.ps_df)<-c("site_branch","fdr_cor_p")
write.table(fdr.ps_df, file=paste0(this_run,".sig-fdr-corrected-ps.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#Perform Bonferroni correction
bon.ps<-p.adjust(split_ps, method="bonferroni")
bon.ps.sig<-bon.ps[which(bon.ps < 0.05)]
print(paste("Length of Bonferroni-corrected P-value vector:", length(bon.ps)))
print(paste("Number of significant P-values (P<0.05) after Bonferroni correction:", length(bon.ps.sig)))
print(head(bon.ps))
bon.ps_df<-as.data.frame(cbind(names(bon.ps.sig), bon.ps.sig))
colnames(bon.ps_df)<-c("site_branch","bon_cor_p")
write.table(bon.ps_df, file=paste0(this_run,".sig-bon-corrected-ps.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#Need to figure out a way to keep these connected to their branch and gene name, then save this information
