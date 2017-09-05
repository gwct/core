#############################################################################
# Wrapper functions -- used by wrappers.py
# Gregg Thomas
# August 2017
#############################################################################

import core, sys, os, subprocess, treeparse as tp

#############################################################################

def ioInfo(input_init, path, output_init, program, file_flag):
# This script defines the output locations for the given options, ensuring that no files
# are overwritten. If a duplicate location is given for the output, and integer count is added
# to that location.
	pad = 40;
	print core.spacedOut("Running " + program + " on FASTA files in:", pad), input_init;
	if path:
		print core.spacedOut("User specified path to " + program + ":", pad), path;
		if not os.path.exists(path):
			sys.exit(core.errorOut(2, "Path (-p) to program not found!!"));
	else:
		print core.spacedOut("No path to " + program + " specified. Trying:", pad), program.lower();
		path = program.lower();
	# Parses the path option.

	output, file_num = defaultOut(input_init, file_flag, program.lower(), output_init);
	# Calls the function that renames the output location based on the previous runs.

	print core.spacedOut("Writing output to:", pad), output;
	if not file_flag:
		print "* Making output directory...";
		os.system("mkdir " + output);
		print "* Making logfile...";
		logfilename = os.path.join(output, "run-" + program.lower() + "-" + str(file_num) + ".log");
	else:
		logfilename = os.path.join(os.path.dirname(output), "run-" + program.lower() + "-" + str(file_num) + ".log");
	core.filePrep(logfilename, "# " + " ".join(sys.argv) + "\n");
	# Makes the output directory (if applicable) and the log file.

	return path, output, logfilename;
	# Returns the parsed path, output location, and log file name.

#############################################################################

def defaultOut(input_name, file_flag, suffix, output_init=False):
# This function decides whether to make a new output location or use
# the user specified one. If the output location already exists it ensures
# that the old location is not overwritten.
	i = 2;
	if file_flag:
	# If the input location is a file.
		if not output_init:
			output = os.path.splitext(input_name);
		# If the user did not specify an output file name, take the base of the input file name.
		else:
			output = os.path.splitext(output_init);
		# Otherwise, use the user specified option.
		output = output[0] + "-" + suffix + "-1" + output[1];

		while os.path.exists(output) or os.path.exists(os.path.splitext(output)[0]):
			output = os.path.splitext(output);
			output = output[0][:output[0].rindex("-")+1] + str(i) + output[1];
			i += 1;
		# If the chosen output file exists, this will continually add 1 to a counter label at the end
		# of the file until a new file is chosen that does not exist.

	else:
	# If the input location is a directory.
		if not output_init:
			output = input_name.rstrip("/").rstrip("\\") + "-" + suffix;
		# If the user did not specify an output name, a directory will be made based on the input directory name.
		else:
			output = output_init;
		# Otherwise, use the user specified option.
		output += "-1";

		while os.path.exists(output):
			output = output[:output.rindex("-")+1] + str(i);
			i += 1;
		# If the chosen output directory exists, this will continually add 1 to a counter label at the end
		# of the directory until a new directory is chosen that does not exist.

	return output, (i-1);

#############################################################################

def runMuscle(infiles, file_flag, path, v, output, logfilename):
# This module runs the MUSCLE alignment program on a list of FASTA files.
	print "Running MUSCLE...\n";
	if v == 0 and not file_flag:
		stdoutlog = os.path.join(output, "muscle.stdout");
	# If the user specifies nothing to be printed to the screen, MUSCLE's output will instead be
	# redirected to a file.

	i, numbars, donepercent, numfiles = 0,0,[], len(infiles);
	fa_skip = [];

	for infile in infiles:
		if v == 0 and not file_flag:
			numbars, donepercent = core.loadingBar(i, numfiles, donepercent, numbars);
		i += 1;
		if not infile.endswith(".fa"):
			fa_skip.append(infile);
			continue;
		# Read the file if it is a FASTA (.fa) file.

		if file_flag:
			outfilename = output;
		else:
			outfilename = os.path.splitext(os.path.basename(infile))
			outfilename = os.path.join(output, outfilename[0] + "-muscle" + outfilename[1]);
		# Get the output file name for the current alignment.

		muscle_cmd = path + " -in '" + infile + "' -out '" + outfilename + "'";
		if v == 0 and not file_flag:
			muscle_cmd += ">> " + stdoutlog + " 2>&1";
		os.system(muscle_cmd);
		if not file_flag:
			with open(logfilename, "a") as logfile:
				logfile.write(muscle_cmd + "\n");
		# The MUSCLE call.

	if v == 0 and not file_flag:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	core.printWrite(logfilename,"-----", file_flag);
	if fa_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(fa_skip)) + " file(s) were skipped because they couldn't be read as fasta files: " + ",".join([os.path.basename(f) for f in fa_skip]), file_flag);
	print "=======================================================================";

#############################################################################

def runPasta(infiles, file_flag, path, seqtype, v, output, logfilename):
# This module runs the PASTA alignment program on a list of FASTA files.
	print "Running PASTA...\n";
	if v == 0 and not file_flag:
		stdoutlog = os.path.join(output, "pasta.stdout");
	# If the user specifies nothing to be printed to the screen, PASTA's output will instead be
	# redirected to a file.

	i, numbars, donepercent, numfiles = 0,0,[], len(infiles);
	fa_skip = [];

	for infile in infiles:
		if v == 0 and not file_flag:
			numbars, donepercent = core.loadingBar(i, numfiles, donepercent, numbars);
		i += 1;
		if not infile.endswith(".fa"):
			fa_skip.append(infile);
			continue;
		# Read the file if it is a FASTA (.fa) file.

		jobname = os.path.splitext(os.path.basename(infile))[0];
		if file_flag:
			pastadir = os.path.splitext(output)[0];
			outfilename = pastadir + ".fa";
		else:
			pastadir = os.path.join(output, jobname + "-pastadir");
			outfilename = pastadir + ".fa";
			outfilename = outfilename.replace("pastadir", "pasta");

		if not os.path.exists(pastadir):
			mkcmd = "mkdir " + pastadir;
			os.system(mkcmd);
		# Get the output info for the current alignment.

		pasta_cmd = "run_pasta.py -i '" + infile + "' -d " + seqtype + " -j " + jobname + " -o '" + pastadir + "'";
		if v == 0 and not file_flag:
			pasta_cmd += ">> " + stdoutlog + " 2>&1";
		os.system(pasta_cmd);
		if not file_flag:
			with open(logfilename, "a") as logfile:
				logfile.write(pasta_cmd + "\n");
		# The PASTA call.

		if file_flag:
			os.system("mv '" + logfilename + "' '" + pastadir + "'");
			logfilename = os.path.join(output, os.path.basename(logfilename));
		else:
			pastafiles = os.listdir(pastadir);
			for pfile in pastafiles:
				if "marker" in pfile and ".aln" in pfile:
					pasta_fasta = os.path.join(pastadir, pfile);
					cpcmd = "cp '" + pasta_fasta + "' '" + outfilename + "'";
					os.system(cpcmd);
					break;
		# If the input type is a file, then keep all the files in the pastadir and move the log file
		# into it as well. Otherwise, move the actual alignment out into the output directory.

	if v == 0 and not file_flag:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	core.printWrite(logfilename,"-----", file_flag);
	if fa_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(fa_skip)) + " file(s) were skipped because they couldn't be read as fasta files: " + ",".join([os.path.basename(f) for f in fa_skip]), file_flag);
	print "=======================================================================";

#############################################################################

def runGblocks(infiles, file_flag, path, seqtype, run_mode, v, output, logfilename):
# This module runs the GBlocks alignment masking program on a list of aligned FASTA files.
	print "Running GBlocks...\n";
	if v == 0 and not file_flag:
		stdoutlog = os.path.join(output, "gblocks.stdout");
	# If the user specifies nothing to be printed to the screen, GBlocks' output will instead be
	# redirected to a file.

	i, numbars, donepercent, numfiles = 0,0,[], len(infiles);
	fa_skip, aln_skip = [],[];
	acc_mask, rej_mask = 0,0;

	for infile in infiles:
		if v == 0 and not file_flag:
			numbars, donepercent = core.loadingBar(i, numfiles, donepercent, numbars);
		i += 1;
		seqs, skip = core.fastaReader(infile);
		if skip:
			fa_skip.append(infile);
			continue;
		if not core.checkAlign(seqs):
			aln_skip.append(infile);
			continue;
		seqlen = len(seqs[seqs.keys()[0]]);
		# Read the file if it is a FASTA (.fa) file and check to make sure it is an alignment.

		gblocks_outfile = infile + "-gb";
		if file_flag:
			outfilename = output;
		else:
			outfilename = os.path.splitext(os.path.basename(infile))
			outfilename = os.path.join(output, outfilename[0] + "-gblocks" + outfilename[1]);
		# Get the output info for the current alignment.

		gblocks_cmd = "gblocks '" + infile + "' -t=" + seqtype;
		if run_mode == 1:
			b1 = int(round(0.5 * len(seqs))) + 1;
			gblocks_cmd += " -b1=" + str(b1) + " -b2=" + str(b1) + " -b3=" + str(seqlen) + " -b4=2 -b5=a";		
		if v == 0 and not file_flag:
			gblocks_cmd += ">> " + stdoutlog + " 2>&1";
		os.system(gblocks_cmd);
		logline = gblocks_cmd;
		# The GBlocks call.

		#####
		percdiff = 100.0;
		gbseqs = core.fastaGetDict(gblocks_outfile);
		gbseqs = {title : gbseqs[title].replace(" ","") for title in gbseqs }
		gblen = len(gbseqs[gbseqs.keys()[0]]);
		difflen = seqlen - gblen;
		percdiff = float(difflen) / float(seqlen) * 100.0;
		logline += "\t" + str(difflen) + " of " + str(seqlen) + " masked (" + str(round(percdiff,2)) + "%)";
		# Calculate the amount of sequence masked relative to the original alignment.

		if percdiff <= 20 or run_mode == 0:
			logline += "\tAccepting masked alignment";
			acc_mask += 1;
			os.system("mv '" + gblocks_outfile + "' '" + outfilename + "'");
		# If run mode is 0, or less than 20% of the sequence was masked, accept the masked alignment and move the
		# GBlocks output file to the output directory.

		else:
			logline += "\tRejecting masked alignment";
			rej_mask += 1;
			if file_flag == 1:
				outfilename = outfilename.replace("gblocks-","gblocks-nomask-");
			else:
				outfilename = outfilename.replace("gblocks.fa","gblocks-nomask.fa");
			os.system("cp '" + infile + "' '" + outfilename + "'");
			os.system("rm '" + gblocks_outfile + "'");
		# If run mode is 1 and more than 20% os the sequence is masked, reject the masked alignment and move the
		# original alignment to the output directory. And remove the GBlocks output files.

		os.system("rm '" + infile + "-gb.htm'");
		# Remove the htm file created by GBlocks. I never use this.

		if not file_flag:
			with open(logfilename, "a") as logfile:
				logfile.write(logline + "\n");
		# My implementation of GBlocks has two run modes: stringent (0) or phylogenetic (1). The stringent mode accepts all
		# masks, while the phylogenetic mode only accepts masks the remove less than 20% of the sequence. This has been shown
		# to be more accurate for phylogenetic reconstructions. If a mask is not accepted, the original alignment is written
		# in its place. This block of code decides whether to Accept or Reject the mask, given the run mode and the amount
		# of sequence masked.	
		#####

	if v == 0:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	core.printWrite(logfilename,"-----", file_flag);
	core.printWrite(logfilename,"# masks accepted:\t" + str(acc_mask), file_flag);
	core.printWrite(logfilename,"# masks rejected:\t" + str(rej_mask), file_flag);
	if fa_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(fa_skip)) + " file(s) were skipped because they couldn't be read as fasta files: " + ",".join([os.path.basename(f) for f in fa_skip]), file_flag);
	if aln_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(aln_skip)) + " file(s) were skipped because they might not have been alignments: " + ",".join([os.path.basename(f) for f in aln_skip]), file_flag);
	print "=======================================================================";

#############################################################################

def runRaxml(infiles, file_flag, path, model, bs_reps, threads, constraint_tree, estimate_bl, v, output, logfilename):
# This module runs RAxML to estimate maximum likelihood phylogenies on a list of FASTA (.fa) alignment files.
	from random import randint
	print "Running RAxML...\n";
	if file_flag:
		output = os.path.splitext(output)[0];
		os.system("mkdir '" + output + "'");
		os.system("mv '" + logfilename + "' '" + output + "'");
		logfilename = os.path.join(output, os.path.basename(logfilename));
	# RAxML creates many output files, so if the input type is a file, we actually want to
	# create an output directory and move everything in there.

	if not file_flag:
		bestdir = os.path.join(output, "raxml-best");
		outdir = os.path.join(output, "raxml-out");
		os.system("mkdir '" + bestdir + "'");
		os.system("mkdir '" + outdir + "'");
	# If we are running RAxML on several files (ie if input type is a directory), then we
	# need to organize the output directory to keep track of everything. 

	if v == 0:
		stdoutlog = os.path.join(output, "raxml.stdout");
	# If the user specifies nothing to be printed to the screen, RAxML's output will instead be
	# redirected to a file.

	seedfile = open(os.path.join(output, "raxml-seeds.txt"), "w");
	if bs_reps > 0:
		bseedfile = open(os.path.join(output, "raxml-bseeds.txt"), "w");
	# RAxML's likelihood search requires a starting seed number. These are chosen randomly and written
	# to a file just to keep track of them.

	i, numbars, donepercent, num_files = 0,0,[],len(infiles);
	fa_skip, no_tree, num_reduced, num_files_read = [],[],0,0;
	trees = {};

	for infile in infiles:
		if v == 0 and not file_flag:
			numbars, donepercent = core.loadingBar(i, num_files, donepercent, numbars);
		i += 1;
		if not infile.endswith(".fa"):
			fa_skip.append(infile);
			continue;
		num_files_read += 1;
		# Read the file if it is a FASTA (.fa) file.

		rax_outfile = os.path.splitext(os.path.basename(infile))[0];
		if not file_flag:
			rax_outdir = os.path.join(outdir, rax_outfile + "-raxout/");
			if not os.path.exists(rax_outdir):
				os.system("mkdir '" + rax_outdir + "'");
		# Get the output info for the current alignment.
		
		seed = str(randint(1000000,999999999));
		seedfile.write(infile + "\t" + str(seed) +"\n");
		if bs_reps > 0:
			boot_seed = str(randint(1000000,999999999));
			bseedfile.write(infile + "\t" + str(boot_seed) +"\n");
		# Generate the starting seed and bootstrap seed (if applicable).

		raxml_cmd = "'" + path + "' ";
		if bs_reps > 0:
			raxml_cmd += "-f a ";
		raxml_cmd += " -m " + model + " -p " + seed;
		if bs_reps > 0:
			raxml_cmd += " -x " + boot_seed + " -# " + str(b);
		if threads > 1:
			raxml_cmd += " -T " + str(t);
		if constraint_tree:
			raxml_cmd += " -g '" + constraint_tree + "'";
			if estimate_bl:
				raxml_cmd += " -f e";
		raxml_cmd += " -s '" + infile + "' -n '" + rax_outfile + "' -w '" + os.path.abspath(output) + "'";
		if v == 0:
			raxml_cmd += " >> " + stdoutlog + " 2>&1";
		# Building the RAxML command based on the input parameters.
		with open(logfilename, "a") as logfile:
			logfile.write(raxml_cmd + "\n");
		os.system(raxml_cmd);
		# The RAxML call

		no_tree_flag = 1;
		if not file_flag:
			raxout_files = os.listdir(output);
			for raxfile in raxout_files:
				full_file = os.path.join(output, raxfile);
				if "RAxML_bestTree" in raxfile:
					if bs_reps == 0:
						trees[rax_outfile] = open(full_file, "r").read();
					os.system("mv '" + full_file + "' '" + bestdir + "'");
					no_tree_flag = 0;
				elif "RAxML" in raxfile:
					if bs_reps > 0 and "bipartitions." in raxfile:
						trees[rax_outfile] = open(full_file, "r").read();
					os.system("mv '" + full_file + "' '" + rax_outdir + "'");
					no_tree_flag = 0;
			if os.path.exists(infile + ".reduced"):
				num_reduced += 1;
				os.system("mv '" + infile + ".reduced' '" + rax_outdir + "'");
		if no_tree_flag == 1:
			no_tree.append(infile);
		# RAxML creates several files for each alignment. We are interested in the "bestTree" file or the
		# "bipartitions" file if bootstrap has been run. These are copied to a separate directory here. The trees
		# are also saved to a dictionary for later. RAxML also creates ".reduced" alignment files that 
		# remove any dubious sites. They are also moved to the output directory for the current run.

	if not file_flag:
		gtfile = open(os.path.join(output, "best-trees.txt"), "w");
		astralfile = open(os.path.join(output, "gt-for-astral.txt"), "w");
		sdmfile = open(os.path.join(output, "gt-for-sdm.txt"), "w");
		if bs_reps > 0:
			astralbsfile = open(os.path.join(output, "bs-for-astral.txt"), "w");

		sdmfile.write(str(len(trees)) + "\n");
		for tree in trees:
			gtfile.write(tree + "\t" + trees[tree]);
			astralfile.write(trees[tree]);
			sdmfile.write(trees[tree]);
			if bs_reps > 0:
				cur_bsfile = os.path.join(rax_outdir, tree + "_raxout", "RAxML_bootstrap." + tree);
				astralbsfile.write(cur_bsfile + "\n");

		gtfile.close();
		astralfile.close();
		sdmfile.close();
		if bs_reps > 0:
			astralbsfile.close();
	# I use the RAxML trees for several downstream analyses (SDM, ASTRAL) and so I write separate files for
	# each.

	if v == 0:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	core.printWrite(logfilename,"-----", file_flag);
	core.printWrite(logfilename,"# of files processed as .fa files:\t" + str(num_files_read), file_flag);
	if fa_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(fa_skip)) + " file(s) were skipped because they couldn't be read as fasta files: " + ",".join([os.path.basename(f) for f in fa_skip]), file_flag);
	core.printWrite(logfilename,"# of files 'reduced' by raxml:\t" + str(num_reduced), file_flag);
	core.printWrite(logfilename,"# of trees made successfully:\t" + str(len(trees)), file_flag);
	if no_tree != []:
		core.printWrite(logfilename,"# The following " + str(len(no_tree)) + " file(s) did not have trees made successfully: " + ",".join([os.path.basename(f) for f in no_tree]), file_flag);
	print "=======================================================================";

#############################################################################

def runCodeml(infiles, file_flag, path, seqtype, treefile, prune, branch_site, anc, v, output, logfilename):
# This module runs codeml for several purposes including the null and alternate models of the branch-site test and for
# ancestral reconstructions of input sequences.
	import treeparse as tr
	print "Running codeml...\n";
	if file_flag:
		output = os.path.splitext(output)[0];
		os.system("mkdir '" + output + "'");
		os.system("mv '" + logfilename + "' '" + output + "'");
		logfilename = os.path.join(output, os.path.basename(logfilename));
	# codeml creates many output files, so if the input type is a file, we actually want to
	# create an output directory and move everything in there.

	if v == 0:
		stdoutlog = os.path.join(output, "codeml.stdout");
	# If the user specifies nothing to be printed to the screen, codeml's output will instead be
	# redirected to a file.

	if anc:
		ancdir = os.path.join(output, "anc-seqs-fa");
		os.system("mkdir '" + ancdir + "'");
	# Create the directory to store ancestral sequences.

	codeml_outdir = os.path.join(output, "codeml-out");
	os.system("mkdir '" + codeml_outdir + "'");
	# Create the codeml-out directory within the main output directory.

	if prune:
		td, tree, root = tr.treeParse(open(treefile, "r").read());
		tips = [node for node in td if td[node][2] == 'tip'];
	# Read the input species tree with treeParse if the --prune option is specified.

	ctlfilename = "codeml.ctl";
	# codeml is executed with a control file.

	i, numbars, donepercent, num_files = 0,0,[],len(infiles);
	fa_skip, aln_skip, num_files_read, num_pruned, no_anc = [],[],0,0,[];

	for infile in infiles:
		if v == 0 and not file_flag:
			numbars, donepercent = core.loadingBar(i, num_files, donepercent, numbars);
		i += 1;
		if not infile.endswith(".fa"):
			fa_skip.append(infile);
			continue;
		# Read the file if it is a FASTA (.fa) file.

		if prune:
			seqs, skip = core.fastaReader(infile);
			if not core.checkAlign(seqs):
				aln_skip.append(infile);
				continue;

			to_prune = [];
			to_prune = [tip for tip in tips if ">" + tip not in seqs];
			if to_prune != []:
				nw_cmd = "nw_prune '" + treefile + "' ";
				for tip in to_prune:
					nw_cmd += tip + " ";
				nw_cmd += " > pruned.tre";
				os.system(nw_cmd);
				num_pruned += 1;
		num_files_read += 1;
		# Use Newick Utilities to prune the input tree to only the species present in the current
		# alignment file.

		codeml_outfile = os.path.splitext(os.path.basename(infile))[0];
		if not file_flag:
			cur_outdir = os.path.join(codeml_outdir, codeml_outfile + "-codemlout/");
			if not os.path.exists(cur_outdir):
				os.system("mkdir '" + cur_outdir + "'");
		# Get the output info for the current alignment.

		with open(ctlfilename, "w") as ctlFile:
			inline = "seqfile = " + infile + "\n";
			ctlFile.write(inline);
			if treefile != "":
				if prune == 1:
					treeline = "treefile = pruned.tre\n";
				else:
					treeline = "treefile = " + treefile + "\n";
				ctlFile.write(treeline);
			if file_flag:
				outline = "outfile = " + os.path.join(output, codeml_outfile + ".out") + "\n\n";
			else:
				outline = "outfile = " + os.path.join(cur_outdir, codeml_outfile + ".out") + "\n\n";
			ctlFile.write(outline);

			ctlFile.write("noisy = 3\n");
			ctlFile.write("verbose = 0\n");
			if treefile != "":
				ctlFile.write("runmode = 0\n\n");
			else:
				ctlFile.write("runmode = 2\n\n");

			if seqtype == 'codon':
				ctlFile.write("seqtype = 1\n");
			elif seqtype == 'protein':
				ctlFile.write("seqtype = 2\n");

			ctlFile.write("CodonFreq = 2\n");
			ctlFile.write("clock = 0\n");
			ctlFile.write("aaDist = 0\n");
			ctlFile.write("aaRatefile = " + os.path.join(path, "dat", "wag.dat") + "\n");
			ctlFile.write("model = 2\n\n");

			if branch_site in [1,2]:
				ctlFile.write("NSsites = 2\n\n");
			else:
				ctlFile.write("NSsites = 0\n\n");

			ctlFile.write("icode = 0\n");
			ctlFile.write("fix_kappa = 0\n");
			ctlFile.write("kappa = 3\n");
			if branch_site == 1:
				ctlFile.write("fix_omega = 1\n");
			else:
				ctlFile.write("fix_omega = 0\n");
			ctlFile.write("omega = 1\n\n");

			ctlFile.write("fix_alpha = 1\n");
			ctlFile.write("alpha = 0\n");
			ctlFile.write("Malpha = 0\n");
			ctlFile.write("ncatG = 10\n\n");

			ctlFile.write("getSE = 0\n");
			if anc:
				ctlFile.write("RateAncestor = 1\n");
			ctlFile.write("Small_Diff = .5e-6\n");
		# Write the control file based on the input options.

		codeml_cmd = os.path.join(path, "bin", "codeml " + ctlfilename);
		if v == 0:
			codeml_cmd += " >> " + stdoutlog + " 2>&1";
		with open(logfilename, "a") as logfile:
			logfile.write(codeml_cmd + "\n");
		os.system(codeml_cmd);
		# The codeml call

		#####
		if anc:
			anc_seqfile = os.path.join(ancdir, codeml_outfile + "-anc.fa");
			anc_probfile = os.path.join(ancdir, codeml_outfile + "-ancprobs.fa");
			anc_treefile = os.path.join(ancdir, codeml_outfile + "-anc.tre");
			# Ancestral sequences are saved to two files and the tree used for this alignment is saved to
			# another.

			rstlines = open("rst", "r").readlines();
			# codeml writes ancestral sequences to the rst output file.

			found_seqs = 0;
			node_list = [];
			for k in xrange(len(rstlines)):
				if rstlines[k] == "tree with node labels for Rod Page's TreeView\n":
					anctree = rstlines[k+1].replace(" ","");
					with open(anc_treefile, "w") as atfile:
						atfile.write(anctree);
				# Find the tree used for this alignment.

				if rstlines[k] == "List of extant and reconstructed sequences\n":
					asfile = open(anc_seqfile, "w");
					found_seqs = 1;
					j = 4;
					while rstlines[k+j].find("Overall accuracy of the") == -1:
						if rstlines[k+j] != "\n":
							tmpline = rstlines[k+j].replace("\n", "");
							tmpline = tmpline.split("  ");
							tmpline = filter(None, tmpline);
							node_list.append(tmpline[0]);

							title = ">" + tmpline[0] + "\n";
							asfile.write(title);
							seq = tmpline[1].replace(" ","") + "\n";
							asfile.write(seq);
						j = j + 1;

					asfile.close();
					break;
				# Find the sequences used and reconstructed for this alignment and write them
				# to the appropriate file in FASTA format.

			curseqs = { n : "" for n in node_list };
			for k in xrange(len(rstlines)):
				if rstlines[k] == "Prob of best state at each node, listed by site\n":
					j = 4;
					while rstlines[k+j] != "\n":
						tmpline = filter(None, rstlines[k+j].replace("\n","").split("  "));
						cur_states = tmpline[2].split(": ");
						extant = list(cur_states[0].replace(" ",""));
						ancs = filter(None, cur_states[1].split(" "));
						all_states = extant + ancs;

						for n in xrange(len(node_list)):
							cur_spec = node_list[n];
							cur_aa = all_states[n];
							if cur_aa.find("(") != -1 and cur_aa.find(")") != -1:
								cur_aa = cur_aa.replace("(","_").replace(")","") + " ";
							curseqs[cur_spec] = curseqs[cur_spec] + cur_aa;
						j = j + 1;

					asfile = open(anc_probfile, "w");
					for seq in curseqs:
						title = ">" + seq + "\n";
						asfile.write(title);
						for base in curseqs[seq]:
							asfile.write(base);
						asfile.write("\n");
					asfile.close();
					break;
			# Find the ancestral sequences with probabilities for this alignment and write them
			# to the appropriate file in FASTA format.

			if found_seqs == 0:
				no_anc.append(infile);
		# If ancestral reconstructions are done, I save them in a specific directory. I save the sequences,
		# the sequences with their probabilities, and the tree with node labels for each alignment.
		#####

	newfilelist = os.listdir(os.getcwd());
	for each in newfilelist:
		if each in ["2NG.dN","2NG.dS","2NG.t","codeml.ctl","lnf","rst","rst1","rub","pruned.tre"]:
			if file_flag:
				mv_cmd = "mv " + each + " " + output;
			else:
				mv_cmd = "mv " + each + " " + cur_outdir;
			os.system(mv_cmd);
	# Move all of codeml's output files to the correct place.

	if v == 0 and not file_flag:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	print "\n" + core.getTime() + " Done!";
	core.printWrite(logfilename,"-----", file_flag);
	core.printWrite(logfilename,"# of files processed as .fa files:\t" + str(num_files_read), file_flag);
	if fa_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(fa_skip)) + " file(s) were skipped because they couldn't be read as fasta files: " + ",".join([os.path.basename(f) for f in fa_skip]), file_flag);
	if aln_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(aln_skip)) + " file(s) were skipped because they might not have been alignments: " + ",".join([os.path.basename(f) for f in aln_skip]), file_flag);
	core.printWrite(logfilename,"# of files that used a pruned species tree:\t" + str(num_pruned), file_flag);
	if anc and no_anc != []:
		core.printWrite(logfilename,"# The following " + str(len(no_anc)) + " file(s) failed to write ancestral sequences: " + ",".join([os.path.basename(f) for f in no_anc]), file_flag);
	print "=======================================================================";

#############################################################################

def runSDM(infiles, file_flag, path, v, output, logfilename):
# This module infers a species tree runs the SDM implementation of average consensus along with the Neighbor-Joining Algorithm as
# implemented in R on a set of gene trees.
	print "Running SDM...\n";
	output = os.path.splitext(output)[0];
	os.system("mkdir '" + output + "'");
	os.system("mv '" + logfilename + "' '" + output + "'");
	logfilename = os.path.join(output, os.path.basename(logfilename));
	# SDM creates many output files, so if the input type is a file, we actually want to
	# create an output directory and move everything in there.

	if v == 0:
		stdoutlog = os.path.join(output, "sdm.stdout");
	# If the user specifies nothing to be printed to the screen, codeml's output will instead be
	# redirected to a file.

	infile = infiles[0];
	indir = os.path.dirname(infile);
	sdmfilename = os.path.join(output, os.path.basename(infile) + "_sdm_mat.txt");
	rmatfilename = os.path.join(output, "rmat.tmp");
	routfilename = os.path.join(output, "sdm-nj.tre");
	# Some more file prep

	sdm_cmd = "java -jar ~/bin/SDM/SDM.jar -i '" + infile + "' -d ACS97 -t T -f Phylip_square";
	if v == 0:
		sdm_cmd += " >> " + stdoutlog + " 2>&1";
	with open(logfilename, "a") as logfile:
		logfile.write(sdm_cmd + "\n");
	os.system(sdm_cmd);
	# The SDM call with ACS97

	os.system("mv '" + os.path.join(indir, infile + "_sdm_deformed_matrices.txt") + "' '" + output + "'");
	os.system("mv '" + os.path.join(indir, infile + "_sdm_mat.txt") + "' '" + output + "'");
	os.system("mv '" + os.path.join(indir, infile + "_sdm_rates.txt") + "' '" + output + "'");
	os.system("mv '" + os.path.join(indir, infile + "_sdm_tab.txt") + "' '" + output + "'");
	# Moves all SDM output files to script outdir

	matlines = open(sdmfilename, "r").readlines();
	specs = [];
	for x in xrange(len(matlines)):
		if x == 0 or matlines[x] == "\n":
			continue;
		while matlines[x].find("  ") != -1:
			matlines[x] = matlines[x].replace("  "," ");
		specs.append(matlines[x][:matlines[x].index(" ")]);

	specs = " ".join(specs);

	with open(rmatfilename, "w") as rmatfile:
		rmatfile.write(specs + "\n");
		for x in xrange(len(matlines)):
			if x == 0 or matlines[x] == "\n":
				continue;
			rmatfile.write(matlines[x]);
	# Re-writes the SDM distance matrix to be readable by R

	r_cmd = "Rscript ~/bin/core/corelib/nj_tree.r '" + rmatfilename + "' '" + routfilename + "'";
	if v == 0:
		r_cmd += " >> " + stdoutlog + " 2>&1";
	with open(logfilename, "a") as logfile:
		logfile.write(r_cmd + "\n");
	os.system(r_cmd);
	##Runs the NJ algorithm within the ape package of R

	if v != 0:
		print "\n ----Unrooted SDM-NJ tree----";
		nj_tree = open(routfilename, "r").read().strip();
		print nj_tree + "\n";
	print "\n" + core.getTime() + " Done!";
	print "=======================================================================";

#############################################################################

def runReights(infiles, file_flag, path, numsites, calspec, calage, output, logfilename):
# This module runs r8s on a species tree to infer branch lengths in units of absolute time.
	print "Running r8s...\n";
	output = os.path.splitext(output)[0];
	os.system("mkdir '" + output + "'");
	os.system("mv '" + logfilename + "' '" + output + "'");
	logfilename = os.path.join(output, os.path.basename(logfilename));
	# We want to save both the r8s input and output files, so if the input type is a file, we actually want to
	# create an output directory and move everything in there.

	infile = infiles[0];
	tree = open(infile, "r").read().strip();
	indir = os.path.dirname(infile);
	rinfilename = os.path.join(output, "r8s-input.txt");
	routfilename = os.path.join(output, "r8s-output.txt");
	# r8s file prep

	calspec = [nodes.split(",") for nodes in calspec.split(" ")];
	calage = calage.split(" ");
	if len(calspec) != len(calage):
		sys.exit(core.errorOut(12, "You must enter the same number of calibration nodes (-calspec) and calibration ages (-calage)"));
	calnodes = [node[0][:3] + node[1][:3] for node in calspec];
	# Parsing of the calibration info.

	with open(rinfilename, "w") as r8sfile:
		r8sfile.write("#NEXUS\nbegin trees;\n");
		outline = "tree in_tree = [&U] " + tree + "\nEnd;\n";
		r8sfile.write(outline);
		r8sfile.write("begin rates;\n");
		outline = "blformat nsites=" + str(numsites) + " lengths=persite ultrametric=no;\n";
		r8sfile.write(outline);
		r8sfile.write("collapse;\n");
		for n in xrange(len(calnodes)):
			outline = "mrca " + calnodes[n] + " " + calspec[n][0] + " " + calspec[n][1] + ";\n";
			r8sfile.write(outline);
		for n in xrange(len(calnodes)):
			if calage[n].find("min") != -1 and calage[n].find("max") != -1:
				curcal = calage[n].split("-");
				outline = "constrain taxon=" + calnodes[n] + " minage=" + curcal[1] + " maxage=" + curcal[3] +";\n";
			elif calage[n].find("min") != -1:
				curcal = calage[n].split("-");
				outline = "constrain taxon=" + calnodes[n] + " minage=" + curcal[1] +";\n";
			elif calage[n].find("max") != -1:
				curcal = calage[n].split("-");
				outline = "constrain taxon=" + calnodes[n] + " minage=" + curcal[1] +";\n";
			else:
				outline = "fixage taxon=" + calnodes[n] + " age=" + calage[n] + ";\n";
			r8sfile.write(outline);
		r8sfile.write("divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;\n");
		r8sfile.write("describe plot=chronogram;\n");
		r8sfile.write("describe plot=tree_description;\n");
		r8sfile.write("end;");
	# Write the r8s input file.

	r8s_cmd = path + " -b -f '" + rinfilename + "' > '" + routfilename + "'";
	os.system(r8s_cmd);
	# The r8s call

	div_tree = open(routfilename, "r").readlines()[-1];
	print "----Smoothed tree----";
	print div_tree;
	print core.getTime() + " Done!";
	print "=======================================================================";

#############################################################################

def runNotung(infiles, file_flag, path, spectree, rearrange, bsthresh, bsroot, v, output, logfilename):
# This module runs Notung on a set of gene trees for bootstrap rearrangement and rooting.
	print "Running Notung...";
	try:
		specline = open(spectree, "r").read();
		td, tree, r = tp.treeParse(specline);
	except:
		sys.exit(core.errorOut(19, "Could not read species tree (-s) as a Newick tree!"));
	# First check to make sure the species tree file is a valid Newick tree.

	if v == 0 and not file_flag:
		stdoutlog = os.path.join(output, "notung.stdout");
	# If the user specifies nothing to be printed to the screen, MUSCLE's output will instead be
	# redirected to a file.

	if file_flag:
		sys.exit();
	# Single file Notung commands not yet implemented

	output = os.path.abspath(output);
	batchfilename = os.path.join("run-notung.batch");
	tre_skip = [];
	with open(batchfilename, "w") as batchfile:
		batchfile.write(spectree + "\n");
		for infile in infiles:
			line = open(infile, "r").read();
			try:
				td, tree, r = tp.treeParse(line);
			except:
				tre_skip.append(infile);
				continue;
			batchfile.write(infile + "\n");
	# This section creates Notungs batch file. The first line contains the file path to the species tree
	# while all subsequent lines contain file paths to gene trees.

	notung_cmd = "java -jar ~/bin/Notung-2.8.1.6-beta/Notung-2.8.1.6-beta.jar -b '" + batchfilename + "'";
	if rearrange:
		notung_cmd += " --rearrange --threshold " + str(bsthresh);
	elif bsroot:
		notung_cmd += " --root";
	notung_cmd += " --speciestag postfix --edgeweights name --treeoutput newick --nolosses --outputdir '" + output + "'";
	if v == 0 and not file_flag:
		notung_cmd += ">> " + stdoutlog + " 2>&1";
	os.system(notung_cmd);
	if not file_flag:
		with open(logfilename, "a") as logfile:
			logfile.write(notung_cmd + "\n");
	# The Notung call

	os.system("mv " + batchfilename + " " + os.path.join(output, batchfilename));
	# Move the batch file to the output directory.

	print "\n" + core.getTime() + " Done!";
	core.printWrite(logfilename,"-----", file_flag);
	core.printWrite(logfilename,"# Total files in input directory:\t" + str(len(infiles)), file_flag);
	if tre_skip != []:
		core.printWrite(logfilename,"# The following " + str(len(tre_skip)) + " file(s) were skipped because they couldn't be read as Newick trees: " + ",".join([os.path.basename(f) for f in tre_skip]), file_flag);
	print "=======================================================================";




