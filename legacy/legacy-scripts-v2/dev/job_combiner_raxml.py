#!/usr/bin/python
import sys, os

if len(sys.argv) not in [2,3] or not os.path.exists(os.path.abspath(sys.argv[1])) or (len(sys.argv) == 3 and sys.argv[2] != "1"):
	print("Usage: $ python job_combiner_2.py [indir] [bootstrap flag = 1]");
	sys.exit();

indir = os.path.abspath(sys.argv[1]);
bs_flag = False;
if len(sys.argv) == 3:
	if sys.argv[2] == "1":
		bs_flag = True;

outdir = os.path.join(indir, "raxml-combined/");
rax_bestdir = os.path.join(outdir, "raxml-best/");
rax_outdir = os.path.join(outdir, "raxml-out/");

if not os.path.isdir(outdir):
	print("+Creating output directory");
	os.system("mkdir " + outdir);
if not os.path.isdir(rax_outdir):
	print("+Creating raxml-out directory");
	os.system("mkdir " + rax_outdir);
if not os.path.isdir(rax_bestdir):
	print("+Creating raxml-best directory");
	os.system("mkdir " + rax_bestdir);

dirlist = os.listdir(indir);
num_files = 0;
tree_fail = [];
trees = {};

print("\n=======================================================================\n");
for jobdir in dirlist:
	if jobdir in [".DS_Store", "raxml-combined", "job_scripts", "job-scripts"] or ".sh" in jobdir:
		continue;
	fulldir = os.path.join(indir, jobdir);
	print("Current job directory: " + fulldir);
	fulllist = os.listdir(fulldir);

	cur_num_files = len([f for f in fulllist if os.path.isfile(os.path.join(fulldir, f))]);
	print("# of files: ", cur_num_files);
	num_files += cur_num_files;

	dir_max = 0;
	raxdir = "";
	for adir in fulllist:
		if "run-raxml" in adir:
			dir_num = int(adir[adir.rindex("-")+1:]);
			if dir_num > dir_max:
				dir_max = dir_num;
				raxdir = adir;
	raxdir = os.path.join(fulldir, raxdir);
	print("Current rax directory: " + raxdir);

	cur_outdir = os.path.join(raxdir, "raxml-out/");
	cur_bestdir = os.path.join(raxdir, "raxml-best/");
	cur_treefile = os.path.join(raxdir, "best-trees.txt");

	print("Getting trees....");
	for line in open(cur_treefile):
		line = line.strip().split("\t");
		geneid = line[0];
		tree = line[1];
		trees[geneid] = tree;

	print("Copying raxml output files...");
	cp_cmd = "cp -r " + cur_outdir + "* " + rax_outdir;
	print(cp_cmd);
	os.system(cp_cmd);

	print("Copying raxml best trees...");
	cp_cmd = cp_cmd = "cp " + cur_bestdir + "* " + rax_bestdir;
	print(cp_cmd);
	os.system(cp_cmd);
		
	print("Checking log...");
	wrapfiles = os.listdir(raxdir);
	logfile = os.path.join(raxdir, [f for f in wrapfiles if '.log' in f][0]);
	for line in open(logfile):
		if "did not have trees made successfully" in line:
			failed = line.split(": ")[1].split(",");
			print("Following files failed:", ",".join(failed));
			tree_fail += failed;

	print("-----------------------------------------");

print("Writing final tree files...");
best_outfile = open(os.path.join(outdir, "best-trees.txt"), "w");
astral_outfile = open(os.path.join(outdir, "gt-for-astral.txt"), "w");
sdm_outfile = open(os.path.join(outdir, "gt-for-sdm.txt"), "w");
sdm_outfile.write(str(len(trees)) + "\n");

if bs_flag:
	astral_bsfile = open(os.path.join(outdir, "bs-for-astral.txt"), "w");

for geneid in trees:
	tree = trees[geneid];
	best_outfile.write(geneid + "\t" + tree + "\n");
	astral_outfile.write(tree + "\n");
	sdm_outfile.write(tree + "\n");

	if bs_flag:
		cur_outdir = os.path.join("raxml-out/", geneid + "-raxout");
		cur_bsfile = os.path.join(cur_outdir, "RAxML_bootstrap." + geneid);
		astral_bsfile.write(cur_bsfile + "\n");

best_outfile.close();
astral_outfile.close();
sdm_outfile.close();
if bs_flag:
	astral_bsfile.close();
print("Done!");
print("-----");
if tree_fail != []:
	print("# The following", len(tree_fail), "file(s) did not have trees made successfully:" + ",".join(tree_fail));
print("=======================================================================\n");