#!/usr/bin/python
import sys, os

# if len(sys.argv) != 6:
# 	print "Usage: $ python job_maker.py [order] [# jobs] [time] [memory] [# ppn]";
# 	sys.exit();

# curorder = sys.argv[1];
# jobs = int(sys.argv[2]);
# time = sys.argv[3];
# mem = sys.argv[4];
# proc = sys.argv[5];

if len(sys.argv) != 10 or "-h" in sys.argv:
	print "Usage: $ python job_maker.py [order] [indir] [outdir] [# jobs] [time] [memory] [# ppn] [# bootstrap] [raxml model]";
	sys.exit();

curorder = sys.argv[1];
indir = os.path.abspath(sys.argv[2]);
outdir = os.path.abspath(sys.argv[3]);
jobs = int(sys.argv[4]);
time = sys.argv[5];
mem = sys.argv[6];
proc = sys.argv[7];
bs = sys.argv[8];
model = sys.argv[9];

if not os.path.isdir(indir):
	print "* ERROR: Specified input directory does not exist.";
	sys.exit();
if not os.path.isdir(outdir):
	print "* MESSAGE: Specified output directory does not exist. Creating it for you.";
	cmd = "mkdir " + outdir;
	print cmd;
	os.system(cmd);

print "Splitting files and writing job scripts...";
filelist = os.listdir(indir);
chunksize = float(len(filelist))/jobs
chunks = [ filelist [i:i + int(chunksize)] for i in range(0, (jobs-1)*int(chunksize), int(chunksize))]
chunks.append(filelist[(jobs-1)*int(chunksize):])

jobdir = os.path.join(outdir, "job-scripts");
if not os.path.isdir(jobdir):
	os.system("mkdir " + jobdir);

jobfiles = [];
chunk_num = 1;
for chunk in chunks:
	print "Chunk", chunk_num;
	chunkdir = os.path.join(outdir, str(chunk_num));
	if not os.path.isdir(chunkdir):
		os.system("mkdir " + chunkdir);
	chunkoutdir = os.path.join(chunkdir, "run-raxml");
	for align in chunk:
		infilename = os.path.join(indir, align);
		outfilename = os.path.join(chunkdir, align);
		os.system("cp " + infilename + " " + outfilename);
	
	jobfilename = os.path.join(jobdir, curorder + "-" + str(chunk_num) + ".pbs");
	jobfiles.append(jobfilename);
	jobfile = open(jobfilename, "w");
	jobfile.write("#!/bin/bash\n");
	jobfile.write("#PBS -k o\n");
	jobfile.write("#PBS -l nodes=1:ppn=" + proc + ",walltime=" + time + ":00:00,vmem=" + mem + "GB\n");
	jobfile.write("#PBS -M grthomas@indiana.edu\n");
	jobfile.write("#PBS -m abe\n");
	jobfile.write("#PBS -N " + curorder + "-" + str(chunk_num) + "\n");
	jobfile.write("#PBS -j oe\n");
	jobfile.write("#PBS -o /N/dc2/scratch/grthomas/i5k/\n");
	jobfile.write("#PBS -d /N/dc2/scratch/grthomas/i5k/\n");

	cmd = "time -p python /N/u/grthomas/Carbonate/bin/core/wrappers.py --raxml -i " + chunkdir + " -p /N/u/grthomas/Carbonate/bin/raxml/raxml-pthreads -model " + model + " ";
	if bs != 0:
		cmd += "-b " + bs + " ";
	if proc != '1':
		cmd += "-t " + str(int(proc)-1) + " ";
	cmd += "-v 0 -o " + chunkoutdir;
	jobfile.write(cmd);

	jobfile.close();
	chunk_num += 1;


submitfilename = os.path.join(outdir, curorder + "-submit.sh");
submitfile = open(submitfilename, "w");
submitfile.write("#!/bin/bash\n");
for jobfile in jobfiles:
	submitfile.write("qsub " + jobfile + "\n");
submitfile.close();

