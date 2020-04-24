#############################################################################
# FASTA functions -- used by fasta_nd_furious.py
# FASTA readers are still located in CORE.
# Gregg Thomas
# August 2017
#############################################################################

import core, sys, os

#############################################################################

def countPos(fasta_files, disp_file=0):
# This function counts the number of sequences and positions in a given set of FASTA files.
	total_files, total_seq, total_pos = 0,0,0;
	fa_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;
		total_files += 1;
		for title in seqs:
			total_seq += 1;
			if disp_file == 1:
				print(title + "\t" + len(seqs[title]));
			total_pos += len(seqs[title]);

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total FASTA files:\t", total_files);
	print("Total sequences:\t", total_seq);
	print("Total positions:\t", total_pos);
	if fa_skip != []:
		print("The following", str(len(fa_skip)), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	print("=======================================================================");

#############################################################################

def countAln(fasta_files, spec_opt):
# This function counts a bunch of stats in a given set of alignments.
	total_aln, total_seq, total_pos, total_col, invariant_sites, variant_sites, sites_with_gaps, sites_all_gaps, total_gaps  = 0,0,0,0,0,0,0,0,0;

	if spec_opt:
		spec_aln = {};

	fa_skip = [];
	aln_skip = [];
	# seqlens = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;
		# Reads the file if it ends with .fa or .fa.gz

		if not core.checkAlign(seqs):
			aln_skip.append(fasta_file);
			continue;
		# Basic check that all the sequences in the file are the same length.

		seqlen = len(seqs[list(seqs.keys())[0]]);
		total_aln += 1;
		total_seq += len(seqs);
		# seqlens.append(str(seqlen));

		# if seqlen < 31:
		# 	print fasta_file, seqlen/3;
		# 	f = os.path.basename(fasta_file);
		# 	cmd = "cp " + fasta_file + " C:/Users/Gregg/Desktop/shortalns/" + f;
		# 	os.system(cmd);

		for i in range(seqlen):
		# For each column in the alignment...
			total_col += 1;
			site = [];
			for title in seqs:
			# ...get all bases in each sequence.
				total_pos += 1;
				base = seqs[title][i];
				if base == "N":
					base = "-";
				if base == "-":
					total_gaps += 1;
				# Converts N to - to be counted as a gap.

				if spec_opt:
				## WILL NEED TO ADD TITLE DELIMITER STUFF.
					if title not in spec_aln:
						spec_aln[title] = [0,0];
						# [total non-gap sites, total gap sites]
					else:
						if base == "-":
							spec_aln[title][1] += 1;
						else:
							spec_aln[title][0] += 1;

				site.append(seqs[title][i]);

			if "-" in site:
				sites_with_gaps += 1;
			site = [base for base in site if base != "-"];
			if site == []:
				sites_all_gaps += 1;
				continue;
			# This filters gaps out of the site before determining if the
			# site is variant or invariant.

			if site.count(site[0]) == len(site):
				invariant_sites += 1;
			else:
				variant_sites += 1;

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total # Seqs\t", total_seq);
	print("Total # Positions\t", total_pos);
	print("Total # Aligns\t", total_aln);
	print("Total # Columns\t", total_col);
	print("# Invariant Sites\t", invariant_sites);
	print("# Variant Sites\t", variant_sites);
	print("Total # Gaps\t", total_gaps);
	print("# Sites with Gaps\t", sites_with_gaps);
	print("# of Sites that are all Gap\t", sites_all_gaps);
	print("-----");
	if spec_opt:
		print("Species counts:");
		for spec in spec_aln:
			print("\t".join([spec, str(spec_aln[spec][0]), str(spec_aln[spec][1])]));
	if fa_skip != []:
		print("The following", str(len(fa_skip)), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	if aln_skip != []:
		print("The following", str(len(aln_skip)), "file(s) were skipped because they might not have been alignments: ", ",".join([os.path.basename(f) for f in aln_skip]));
	print("=======================================================================");

	# with open("alnlens.txt", "w") as outfile:
	# 	outfile.write(",".join(seqlens));
#############################################################################

def concat(fasta_files, header_delim, outfilename):
# This function takes a list of FASTA files and concatenates the sequences with common headers.
	partfilename = os.path.splitext(outfilename)[0] + "-partitions.txt";
	# Also writes a partitions file to keep track of which positions came from which files.

	specs = [];
	concats = {};
	# specs will contain just the list of common headers while concats will be [header]:[concatenated sequence].

	total_aln, total_seq = 0,0;
	fa_skip = [];
	aln_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;
		if not core.checkAlign(seqs):
			aln_skip.append(fasta_file);
			continue;

		seqlen = len(seqs[list(seqs.keys())[0]]);
		total_aln += 1;
		for title in seqs:
			total_seq += 1;
			if header_delim != False:
				if header_delim not in title:
					print("\n** WARNING! Specified header delimiter (-delim) not found for the following sequence:\t" + os.path.basename(fasta_file) + " : " + title);
					sys.exit("Exiting...\n");
				else:
					cur_spec = title[1:title.index(header_delim)];
			else:
				cur_spec = title[1:];
			if cur_spec not in specs:
				specs.append(cur_spec);
			if cur_spec not in concats:
				concats[">" + cur_spec] = "";
	# This first loop has to read through all alignments to get all the common headers. In alignments without a
	# given header, N's must be added for that header.

	for each in aln_skip+fa_skip:
		if each in fasta_files:
			fasta_files.remove(each);
	# Remove any files from the list that couldn't be read as FASTA files.

	pfile = open(partfilename, "w");
	pfile.write("");
	total_pos = 0;
	# Open the partitions file.

	#i = 0;
	#numbars = 0;
	#donepercent = [];

	for fasta_file in fasta_files:
		#numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
		#i = i + 1;
		seqs, skip = core.fastaReader(fasta_file);
		seqlen = len(seqs[list(seqs.keys())[0]]);
		pfile.write(fasta_file + " = " + str(total_pos+1) + "-" + str(total_pos+seqlen) + "\n");
		total_pos = total_pos + seqlen;
		# Read the sequences and add the file name and positions to the partitions file.

		cur_specs = [];
		for title in seqs:
			if header_delim != False:
				cur_spec = title[1:title.index(header_delim)];
			else:
				cur_spec = title[1:];
			cur_specs.append(cur_spec);
			concats[">" + cur_spec] += seqs[title];
		# Gets the sequences for all headers in the current alignment and adds them to
		# the concatenated alignment.

		for spec in specs:
			if spec not in cur_specs:
				concats[">" + spec] += "N" * seqlen;
		# If any headers weren't found in the current alignment, add N's for the length
		# of the current alignment for that header.

	pfile.close();
	outfile = open(outfilename, "w");
	for spec in concats:
		outfile.write(spec + "\n");
		outfile.write(concats[spec] + "\n");
	outfile.close();
	# Write the concatenated alignments to the output file.

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total alignments:\t", total_aln);
	print("Total sequences:\t", total_seq);
	print("Total columns:\t", total_pos);
	if fa_skip != []:
		print("The following", str(len(fa_skip)), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	if aln_skip != []:
		print("The following", str(len(aln_skip)), "file(s) were skipped because they might not have been alignments: ", ",".join([os.path.basename(f) for f in aln_skip]));
	print("=======================================================================");
	
#############################################################################

def combine(fasta_files, outfilename):
# This function just takes a bunch of sequences from several FASTA files and puts them
# into a single FASTA file.
	total_files, total_seq, total_pos = 0,0,0;
	fa_skip = [];
	outfile = open(outfilename, "w");
	for fasta_file in fasta_files:
		print("Loading file", fasta_file + "...");
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;
		total_files += 1;
		for title in seqs:
			total_seq += 1;
			total_pos += len(seqs[title]);
			outfile.write(title + "\n");
			outfile.write(seqs[title] + "\n");
	outfile.close();
	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total FASTA files:\t", total_files);
	print("Total sequences:\t", total_seq);
	print("Total positions:\t", total_pos);
	if fa_skip != []:
		print("The following", str(len(fa_skip)), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	print("=======================================================================");

#############################################################################

def split(fasta_files, header_delim, outdir):
# This function takes a single FASTA file with multiple sequences and splits them
# into individual files.
	fasta_file = fasta_files[0];
	if not os.path.isdir(outdir):
		print("\n++ Creating output directory...");
		os.system("mkdir " + outdir);

	total_seq = 0;
	seqs, skip = core.fastaReader(fasta_file);
	if skip:
		sys.exit(core.errorOut(6, "Something went wrong when reading your input file! Does it have the .fa extension? Is it a properly formatted FASTA file?"))

	for title in seqs:
		total_seq += 1;
		if not header_delim:
			outfilename = os.path.join(outdir, title[1:] + ".fa");
		else:
			outfilename = os.path.join(outdir, title[1:title.index(header_delim)] + ".fa");
		outfile = open(outfilename , "w");
		outfile.write(title + "\n");
		outfile.write(seqs[title]);
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");	
	print("Total sequences:\t", total_seq);
	print("=======================================================================");

#############################################################################

def trim(fasta_files, header_delim, file_flag, out_dest):
# This function trims the FASTA headers at the first occurrence of a given string.
	fa_skip = [];
	header_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;

		outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "trim");
		outfile = open(outfilename, "w");
		for title in seqs:
			new_title = title;
			if header_delim not in title:
				if fasta_file not in header_skip:
					header_skip.append(fasta_file);
			else:
				new_title = title[:title.index(header_delim)];
			outfile.write(new_title + "\n");
			outfile.write(seqs[title] + "\n");
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	if fa_skip != []:
		print("The following", str(len(fa_skip)), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	if header_skip != []:
		print("The following", str(len(header_skip)), "file(s) had headers without the delimeter. These sequences were written as is: ", ",".join([os.path.basename(f) for f in header_skip]));
	print("=======================================================================");

#############################################################################

def relabel(fasta_files, ropt, header_delim, new_labels, file_flag, out_dest):
# This function either adds a label to the beginning or end of a FASTA header or 
# completely replaces the FASTA header.
	if header_delim == False:
		header_delim = "_";

	new_labels = new_labels.split(",");
	if len(new_labels) == 1:
		sep_labels = False;
		new_label = new_labels[0];
	else:
		sep_labels = True;
		label_dict = {};
		try:
			for label in new_labels:
				label = label.split(":");
				label_dict[label[0]] = label[1];
		except:
			sys.exit(core.errorOut(9, "There is something wrong with your new label format (-newlabel)!"));
	# Parses the label format;

	total_files = 0;
	fa_skip = [];
	label_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;

		total_files += 1;
		outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "relab");
		outfile = open(outfilename, "w");
		for title in seqs:
			if not sep_labels:
				new_title = core.relabelHeader(title, new_label, header_delim, ropt);
				outfile.write(new_title + "\n");

			elif sep_labels:
				label_found = False;
				for old_label in label_dict:
					if old_label in title:
						label_found = True;
						new_label = label_dict[old_label];
						new_title = core.relabelHeader(title, new_label, header_delim, ropt);
						outfile.write(new_title + "\n");
						break;
				if label_found == False:
					label_skip.append((fasta_file, title));
			outfile.write(seqs[title] + "\n");
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total files written:\t", total_files)
	if fa_skip != []:
		print("The following", len(fa_skip), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	if label_skip != []:
		print("The following file(s) had headers that didn't contain any old labels. They were re-written as they are in the new file.")
		for each in label_skip:
			print(os.path.basename(each[0]) + "\t" + each[1]);
	print("=======================================================================");

#############################################################################

def removeSeq(fasta_files, labels, file_flag, out_dest):
# This function scans the headers of all FASTA files given and removes 
# sequences with a certain label.
	labels = labels.split(",");
	total_files, total_seq, total_seq_rm = 0,0,0;
	fa_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;

		total_files += 1;
		outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "rmseq");
		outfile = open(outfilename, "w");
		for title in seqs:
			total_seq += 1;
			if any(label in title for label in labels):
				total_seq_rm += 1;
				continue;
			outfile.write(title + "\n");
			outfile.write(seqs[title] + "\n");
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total files read:\t", total_files);
	print("Total sequences read:\t", total_seq);
	print("Total sequences removed:\t", total_seq_rm);
	if fa_skip != []:
		print("The following", len(fa_skip), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	print("=======================================================================");

#############################################################################

def removeStarts(fasta_files, seqtype, file_flag, out_dest):
# This function removes either start codons or start M's, depending on the
# type of input sequence specified.
	total_files, total_seq, total_start_rm = 0,0,0;
	fa_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;

		total_files += 1;
		outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "rmstart");
		outfile = open(outfilename, "w");
		for title in seqs:
			seq = seqs[title];
			total_seq += 1;
			if seqtype[0] == 'p' and seq[0].upper() == "M":
				total_start_rm += 1;
				seq = seq[1:];
			elif seqtype[0] == 'c' and seq[:3].upper() == "ATG":
				total_start_rm += 1;
				seq = seq[3:];
			outfile.write(title + "\n");
			outfile.write(seq + "\n");
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total files read:\t", total_files);
	print("Total sequences read:\t", total_seq);
	print("Total starts removed:\t", total_start_rm);
	if fa_skip != []:
		print("The following", len(fa_skip), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	print("=======================================================================");

#############################################################################

def removeStops(fasta_files, seqtype, file_flag, out_dest):
# This function removes either start codons or start M's, depending on the
# type of input sequence specified.
	stop_codons = ["TAG", "TAA", "TGA", "UAG", "UAA", "UGA"];
	total_files, total_seq, total_stop_rm = 0,0,0;
	fa_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;

		total_files += 1;
		outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "rmstop");
		outfile = open(outfilename, "w");
		for title in seqs:
			seq = seqs[title];
			total_seq += 1;
			if seq[-3:] in stop_codons:
				seq = seq[:-3];
				total_stop_rm += 1;
			outfile.write(title + "\n");
			outfile.write(seq + "\n");
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total files read:\t", total_files);
	print("Total sequences read:\t", total_seq);
	print("Total starts removed:\t", total_stop_rm);
	if fa_skip != []:
		print("The following", len(fa_skip), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	print("=======================================================================");

#############################################################################

def replaceBase(fasta_files, replacements, file_flag, out_dest):
# This function reads all sequences in the input and replaces a given state with
# another given state.
	replacements = replacements.split(",");
	replacements = [r.upper() for r in replacements];

	total_files, total_seq, total_pos, total_repl = 0,0,0,0;
	fa_skip = [];
	for fasta_file in fasta_files:
		seqs, skip = core.fastaReader(fasta_file);
		if skip:
			fa_skip.append(fasta_file);
			continue;

		total_files += 1;
		if any("*" in r for r in replacements) or any(" " in r for r in replacements):
			outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "repl");
		else:
			outfilename = core.getOutFile(fasta_file, file_flag, out_dest, "repl." + ".".join(replacements));
		outfile = open(outfilename, "w");
		for title in seqs:
			seq = seqs[title].upper();
			total_seq += 1;
			total_pos += len(seq);
			for replacement in replacements:
				total_repl += seq.count(replacement[0]);
				seq = seq.replace(replacement[0],replacement[1]);
			outfile.write(title + "\n");
			outfile.write(seq + "\n");
		outfile.close();

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total files read:\t", total_files);
	print("Total sequences read:\t", total_seq);
	print("Total positions read:\t", total_pos);
	print("Total replacements made:\t", total_repl);
	if fa_skip != []:
		print("The following", len(fa_skip), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	print("=======================================================================");

#############################################################################

def extractSeqs(fasta_files, titles, delim, file_flag, out_dest):
# This sequence cycles through all input files and extracts sequences if
# they are in the input title list.

	with open(out_dest, "w") as outfile:
		total_files, total_seq = 0,0;
		fa_skip, ext = [], [];
		for fasta_file in fasta_files:
			seqs, skip = core.fastaReader(fasta_file, "ind");
			if skip:
				fa_skip.append(fasta_file);
				continue;

			total_files += 1;
			for title in seqs:
				if delim:
					t = title[1:title.index(delim)];
				else:
					t = title[1:];

				if t in titles:
					seq = core.fastaGetInd(fasta_file, seqs[title])[1];
					ext.append(t);
					outfile.write(">" + t + "\n");
					outfile.write(seq + "\n");

	print("\n" + core.getTime() + " Done!");
	print("-----");
	print("Total files read:\t", total_files);
	print("Total sequences read:\t", total_seq);
	print("Total sequences extracted:\t", len(ext));
	if fa_skip != []:
		print("The following", len(fa_skip), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	if set(titles) != set(ext):
		print("The following sequences were not extracted:\t" + ",".join([t for t in titles if t not in ext]));
	print("=======================================================================");					
