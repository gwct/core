#!/usr/bin/python
########################################################################################
#Script to filter out isoforms from peptide files in FASTA ENSEMBL or NCBI format. This
#script can also add an easier to read species label to each sequence within the file.
#
#Sample ENSEMBL usage: python isoform_filter.py -i [input_fasta_file] -t ens -l [species_label] -o [output_filename]
#
#Sample NCBI usage: python isoform_filter.py -i [input_fasta_file] -t ncbi -g [toplevel_gff_file] -l [species_label] -o [output_filename]
#
#To just do the relabeling, set -f 0. You shouldn't need a gff file for the NCBI file in
#this case. For NCBI relabeling, the gene ID is also moved to the front of the title line.
#
#Written by: Gregg Thomas, Summer 2014
#
#NCBI filter command kindly provided by the folks at NCBI.
#
########################################################################################

import sys, re, os, argparse
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/")
import core

############################################
#Function definitions.

def optParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser()

	parser.add_argument("-i", dest="input_file", help="An input file containing peptides from a species in FASTA format");
	parser.add_argument("-t", dest="file_type", help="Currently supported file types are ENSEMBL and NCBI peptide files. Enter as 'ens' or 'ncbi' here. Note: If file type is NCBI you will also need to specify the top level gff file with -g")
	parser.add_argument("-g", dest="gff_file", help="If file type is NCBI, the top level gff file is also needed and should be specified here.");
	parser.add_argument("-l", dest="spec_label", help="A species label to add to the gene ID of each sequence.", default="");
	parser.add_argument("--cds", dest="cds_opt", help="Set if input file is a CDS file from NCBI.", action="store_true", default=False);
	parser.add_argument("-o", dest="output_file", help="The desired name of the output file. If none is specified the default is [input_filename]_isofiltered.fa or [input_filename]_isofiltered_relabel.fa");

	args = parser.parse_args();

	args.file_type = args.file_type.lower();

	if None in [args.input_file, args.file_type, args.output_file]:
		sys.exit(core.errorOut(1, "An input file (-i), input file type (-t), and an output file (-o) must all be specified."));

	if args.file_type not in ['ens', 'ncbi']:
		sys.exit(core.errorOut(2, "File type (-t) must be one of either 'ens' (Ensembl) or 'ncbi'."));

	if args.file_type == "ens" and args.gff_file != None:
		sys.exit(core.errorOut(3, "A gff file (-g) should not be specified with file type ens."));

	if args.file_type == "ncbi" and args.gff_file == None:
		sys.exit(core.errorOut(4, "A gff file (-g) must be specified with file type ncbi."));

	if args.cds_opt and args.file_type != "ncbi":
		sys.exit(core.errorOut(5, "CDS files (--cds) only accepted with NCBI input (-t ncbi)."));

	if args.spec_label != "":
		args.spec_label += "_";

	return args.input_file, args.file_type, args.gff_file, args.spec_label, args.cds_opt, args.output_file;

############################################
def ensFilter(inseqs, spec_label, outfilename):

	print("Indexing", len(inseqs), "sequences to be filtered.");
	print("Parsing identifiers...");

	for title in inseqs:
		geneid = title[title.index("gene:") + 5:title.index("gene:") + 23];

		if geneid in identDict:
			identDict[geneid].append((title, inseqs[title]));
		else:
			identDict[geneid] = [];
			identDict[geneid].append((title, inseqs[title]));
	sys.stderr.write('\b');

	print("Filtering and writing sequences...");
	numbars, donepercent, i = 0,[],0;

	for key in identDict:
		numbars, donepercent = core.loadingBar(i, len(identDict), donepercent, numbars);

		if len(identDict[key]) == 1:
			long_title, long_seq = identDict[key][0];

		else:
			titlelist = [];
			seqlist = [];

			for tup in identDict[key]:
				cur_itle, cur_seq = tup;
				titlelist.append(cur_itle);
				seqlist.append(cur_seq);

			long_seq = max(seqlist, key=len)
			long_title = titlelist[seqlist.index(long_seq)];

		new_title = ">" + spec_label + long_title[1:];
		core.writeSeq(outfilename, long_seq, new_title);
		i += 1;

	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
	print("\nDone!");
	print(i, "sequences written.");
	print(len(inseqs) - i, "sequences filtered.");

############################################
def ncbiFilter(inseqs, gff_file, spec_label, cds_opt, outfilename):

	numbars, donepercent, i = 0, [], 0;


	print("Obtaining longest isoforms from .gff file...");

	cmd = "zcat " + gff_file + " | awk \'BEGIN{FS=\"	\";OFS=\"|\"}$3==\"CDS\"{if($4<$5){print $5-$4+1,$9}else{print $4-$5+1,$9}}\' | grep \"[NX]P[_]\" | sed \'s/\([0-9]*\).*GeneID:\([0-9]*\).*\([NX]P[_][0-9]*\.[0-9]*\).*/\\1|\\2|\\3/\' | awk \'BEGIN{FS=\"|\";OFS=\"\t\";gene=\"\";acc=\"\";len=0}{if(acc!=$3){print gene,acc,len/3-1;gene=$2;acc=$3;len=$1}else{len=len+$1}}END{print gene,acc,len/3-1}\' | sort -k1,1n -k3,3nr -k2,2 | awk \'BEGIN{FS=\"	\";OFS=\"	\";gene=\"\";acc=\"\";len=0}{if(gene!=$1){print $1,$2,$3};gene=$1;acc=$2;len=$3}\' > ncbi_isoform_filter_tmp11567.txt"
	os.system(cmd);

	tmpFile = open("ncbi_isoform_filter_tmp11567.txt", "r");
	tmpLines = tmpFile.readlines();
	tmpFile.close();
	os.system("rm ncbi_isoform_filter_tmp11567.txt");

	longest_isos = [];

	for each in tmpLines:
		longest_isos.append(each.split("\t")[1]);
	longest_isos = [_f for _f in longest_isos if _f];

	print("Writing longest isoforms to output file...");

	count = 0;

	for title in inseqs:
		numbars, donepercent = core.loadingBar(i, len(inseqs), donepercent, numbars);
		i += 1;

		found = 0;

		for gid in longest_isos:
			if gid in title:
				if cds_opt or "|" not in title:
					new_title = ">" + spec_label + title[1:];
				else:
					gid = title[title.index("P_")-1:title.index("|",title.index("P_"))]
					new_title = ">" + spec_label + gid + " |" + title[1:title.index("P_")-1] + title[title.index("|",title.index("P_"))+1:];	
				core.writeSeq(outfilename, inseqs[title], new_title);
				count += 1;
				break;

	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
	print("\nDone!");
	print(count, "sequences written.");
	print(len(inseqs) - count, "sequences filtered.");

############################################
#Main Block
############################################

infilename, in_type, gff_file, label, cds_opt, outfilename = optParse();

pad = 50;
print("=======================================================================");
print("\t\t\t" + core.getDateTime());
print(core.spacedOut("Filtering isoforms from:", pad) + infilename);
if in_type == "ens":
	print(core.spacedOut("File type:", pad) + "Ensembl");
if in_type == "ncbi":
	print(core.spacedOut("File type:", pad) + "NCBI");
	print(core.spacedOut("Using GFF file:", pad) + gff_file);
	if cds_opt:
		print(" + CDS input");
if in_type == "crow":
	print(core.spacedOut("File type:", pad) + "Crow");
if label != "":
	print(core.spacedOut("Adding label to beginning of FASTA headers:", pad) + label);
print(core.spacedOut("Writing output to:", pad) + outfilename);
core.filePrep(outfilename);
print("--------------------------");

identDict = {};
ins, skip_flag = core.fastaReader(infilename);

if in_type == "ens":
	ensFilter(ins, label, outfilename);
elif in_type == "ncbi":
	ncbiFilter(ins, gff_file, label, cds_opt, outfilename);

print("=======================================================================");









## DEFUNCT FILTER FOR THE CROW FILES
# elif in_type == "crow":
# 	crowFilter(ins, label, outfilename);

# def crowFilter(inSeqs, filterflag, speclabel, outFilename):
# 	rotator = 0;
# 	numbars = 0;
# 	donepercent = [];
# 	i = 0;

# 	if filterflag == 1:
# 		print "Indexing", len(inSeqs), "sequences to be filtered.";
# 		print "Parsing identifiers...";

# 		for each in inSeqs:

# 			rotator = core.loadingRotator(i, rotator, 100)

# 			curTitle, curSeq = core.getFastafromInd(inFilename, each[0], each[1], each[2], each[3]);

# 			if "gene=" not in curTitle:
# 				print curTitle;
# 				continue;

# 			geneid = curTitle[curTitle.index("gene=") + 5:].strip();

# 			if geneid in identDict:
# 				identDict[geneid].append(each);

# 			else:
# 				identDict[geneid] = [];
# 				identDict[geneid].append(each);

# 			i = i + 1;

# 		sys.stderr.write('\b');

# 		print "Filtering and writing sequences...";

# 		i = 0;
# 		count = 0;

# 		for key in identDict:

# 			numbars, donepercent = core.loadingBar(i, len(identDict), donepercent, numbars);

# 			if len(identDict[key]) == 1:
# 				curTitle, curSeq = core.getFastafromInd(inFilename, identDict[key][0][0], identDict[key][0][1], identDict[key][0][2], identDict[key][0][3]);

# 				if speclabel != "":
# 					newTitle = ">" + speclabel + "_" + curTitle[1:];
# 					core.writeSeq(outFilename, curSeq, newTitle);
# 				else:
# 					core.writeSeq(outFilename, curSeq, curTitle);

# 				count = count + 1;

# 			else:
# 				titlelist = [];
# 				seqlist = [];

# 				for inds in identDict[key]:
# 					aTitle, aSeq = core.getFastafromInd(inFilename, inds[0], inds[1], inds[2], inds[3]);

# 					titlelist.append(aTitle);
# 					seqlist.append(aSeq);

# 				longseq = max(seqlist, key=len)

# 				for inds in identDict[key]:
# 					aTitle, aSeq = core.getFastafromInd(inFilename, inds[0], inds[1], inds[2], inds[3]);

# 					if aSeq == longseq:
# 						curTitle, curSeq = core.getFastafromInd(inFilename, inds[0], inds[1], inds[2], inds[3]);

# 						if speclabel != "":
# 							newTitle = ">" + speclabel + "_" + curTitle[1:];
# 							core.writeSeq(outFilename, curSeq, newTitle);
# 						else:
# 							core.writeSeq(outFilename, curSeq, curTitle);

# 						count = count + 1;
# 						break;

# 			i = i + 1;

# 		pstring = "100.0% complete.";
# 		sys.stderr.write('\b' * len(pstring) + pstring);
# 		print "\nDone!";
# 		print count, "out of", len(identDict), "identifiers written.";
# 		print len(inSeqs) - count, "sequences filtered.";

# 	else:
# 		print "Relabeling...";
# 		for seq in inSeqs:
	
# 			numbars, donepercent = core.loadingBar(i, len(inSeqs), donepercent, numbars);
# 			i = i + 1;

# 			curTitle, curSeq = core.getFastafromInd(inFilename, seq[0], seq[1], seq[2], seq[3]);

# 			newTitle = ">" + speclabel + "_" + curTitle[1:];

# 			core.writeSeq(outFilename, curSeq, newTitle);


# 		pstring = "100.0% complete.";
# 		sys.stderr.write('\b' * len(pstring) + pstring);
# 		print "\nDone!";

