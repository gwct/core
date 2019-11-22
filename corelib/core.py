#############################################################################
#CORE functions
#Gregg Thomas
#August 2013-present
#############################################################################

import string, sys, os, re, subprocess, datetime, gzip
from collections import defaultdict

#############################################################################

def dnaCheck(seq, i_name):
#dnaCheck does several things to ensure the input sequence is in proper format (DNA FASTA).

	dna_symbols = ["A","T","C","G","N","-"];

	seq = seq.upper();
	#It first converts the sequence to all uppercase letters.
	seq = seq.replace('U', 'T');
	#The function then converts any U's to T's, changing RNA to DNA.

	for n in range(len(seq)):
	#This loop then reads through every letter of the sequence and makes sure they are all nucleotide symbols.
	#If they are not, an error is printed and the program exits.
		if seq[n] not in dna_symbols:		
			print("\nError! Input sequence from the file", i_name, "is not a nucleotide sequence. Please input only DNA sequences.\n")
			sys.exit();
	return seq;
	#If the input sequence checks out, it is returned to the call.

#############################################################################

def loadingBar(counter, length, done, bars, firstbar=False, disperc=False):
#This function serves as a text loading bar for long scripts with counters. The following
#lines must be added within the script to initialize and terminate the script:
#Initilization:
#numlines = core.getFileLen(alnfilename);
#numbars = 0;
#donepercent = [];
#i = 0;
#Termination:
#	pstring = "100.0% complete.";
#	sys.stderr.write('\b' * len(pstring) + pstring);
#	print "\nDone!";
#
#If length is lines in a file use the core.getFileLen function to get the number of lines in the file

	# try:
	# 	if sys.version[0] == '2':
	# 		pchr = u'\u2591'.encode('utf-8');
	# 		lchr = u'\u2588'.encode('utf-8');
	# 	elif sys.version[0] == '3':
	# 		pchr = u'\u2591';
	# 		lchr = u'\u2588';		
	# except:
	# 	pchr, lchr = "*","*";

	try:
		#pchr, lchr=u'\u2591',u'\u2588';
		pchr, lchr=u'\u2588',u'\u2591';
	except:
		pchr, lchr = "=","|";

	percent = float(counter) / float(length) * 100.0;
	percentdone = int(percent);

	p = str(percent)
	pstring = " " + p[:5] + "% complete.";

	if percentdone % 2 == 0 and done != None and percentdone not in done:
		loading = "|";
		j = 0;
		while j < bars:
			loading += pchr;
			j += 1;
		if j <= 49:
			loading += lchr;
		else:
			loading += pchr;
		j += 1;
		if j == 50:
			loading = loading[:-1] + pchr;

		while j < 50:
			loading += "-";
			j += 1;
		loading += "|";

		if disperc:
			loading += "                 ";
		if firstbar:
			sys.stderr.write(loading);
			firstbar = False
		else:
			sys.stderr.write('\b' * len(loading) + loading);

		done.append(percentdone);
		bars = bars + 1;
	if disperc:
		sys.stderr.write('\b' * len(pstring) + pstring);
	sys.stderr.flush();
	
	return bars, done, firstbar;

#############################################################################

def loadingRotator(counter, rotate, divisor):
#Provides a loading rotator for loops. The following line must be used to initialize the function
#before the loop in the main code:
#rotator = 0;

	rotation = ['|', '/', '-', '\\'];

	if counter % divisor == 0:
		sys.stderr.write('\b' + rotation[rotate]);
		rotate = rotate + 1;
		if rotate >= len(rotation):
			rotate = 0;

	return rotate;

#############################################################################

def bioTranslator(i_name):
#This function takes a FILE of DNA sequences in FASTA format and translates it
#to the corresponding AA sequence. It then returns that sequence as a dictionary.
#
#This function may become deprecated as the call to fastaGetDict can be performed
#before the call to this function, increasing functionality. See newbioTranslator
#below.

	codons = [''] * 5;
	codons[0] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	codons[1] = "---M---------------M---------------M----------------------------";
	codons[2] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	codons[3] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	codons[4] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

	dnaSeqs = fastaGetDict(i_name);
	aaSeqs = {};

	for seq in dnaSeqs:
		aaSeqs[seq] = "";
		if len(dnaSeqs[seq]) % 3 == 0:
			end = len(dnaSeqs[seq]) - 3;		
			k = 0;
			while k <= end:
				currentCodon = dnaSeqs[seq][k:k+3]

				if len(currentCodon) < 3:
					currentAA = ' ';
				elif "-" in currentCodon or "N" in currentCodon or "?" in currentCodon:
					currentAA = "X";
				elif "*" in currentCodon:
					currentAA = "*";
				else:
					first = codons[2].index(currentCodon[0]);
					second = codons[3].index(currentCodon[1], first);
					final = codons[4].index(currentCodon[2], second);

					currentAA = codons[0][final];

				aaSeqs[seq] = aaSeqs[seq] + currentAA;
				k = k + 3;

	return aaSeqs

#x[y] = re.sub(r'\s', '', x[y]);
#This is a simple line to remove all return characters from a string. Must remember to import re above.
#############################################################################

def newbioTranslator(seq):
#This function takes a DNA sequence as a single string and returns the
#corresponding AA sequence.

	codons = [''] * 5;
	codons[0] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	codons[1] = "---M---------------M---------------M----------------------------";
	codons[2] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	codons[3] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	codons[4] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

	aaSeq = "";

	if len(seq) % 3 == 0:
		end = len(seq) - 3;		
		k = 0;

		while k <= end:
			currentCodon = seq[k:k+3]

			if len(currentCodon) < 3:
				currentAA = ' ';
			elif "-" in currentCodon or "N" in currentCodon or "?" in currentCodon:
				currentAA = "X";
			elif "*" in currentCodon:
				currentAA = "*";
			else:
				first = codons[2].index(currentCodon[0]);
				second = codons[3].index(currentCodon[1], first);
				final = codons[4].index(currentCodon[2], second);

				currentAA = codons[0][final];

			aaSeq = aaSeq + currentAA;
			k = k + 3;

	return aaSeq;

#############################################################################

def variance(data):
#Calculates and returns the variance of a list of numbers.
	mean = 0.0;
	for d in data:
		mean = mean + float(d);
	mean = mean / float(len(data));
	var = 0.0;
	for d in data:
		var = var + (float(d) - mean)**2;
	var = var / ((float(len(data)) - 1.0));
	return var;

#############################################################################

def getFileLen(i_name):
#Calls 'wc -l' to get the number of lines in the file.
	p = subprocess.Popen(['wc', '-l', i_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	result, err = p.communicate();
	if p.returncode != 0:
		raise IOError(err);
	return int(result.strip().split()[0]);

#############################################################################

def getDateTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

#############################################################################

def getTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%I:%M:%S");

#############################################################################

def getLogTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y-%I.%M.%S");

#############################################################################

def printWrite(o_name, o_line, file_flag=True):
# Function to print a string AND write it to the file.
	print(o_line);
	if file_flag == False:
		with open(o_name, "a") as f:
			f.write(o_line + "\n");

#############################################################################

def PW(o_line, o_name, file_flag=True):
# Function to print a string AND write it to the file.
	print(o_line);
	if file_flag == True:
		with open(o_name, "a", encoding="utf-8") as f:
			f.write(o_line + "\n");

#############################################################################

def logCheck(lopt, lfilename, outline):
# Function checks whether or not to write to a logfile, print something, or both.
	if lopt == 1:
		printWrite(lfilename, outline);
	else:
		print(outline);

#############################################################################

def errorOut(errnum, errmsg, ropt=0):
# Formatting for error messages.
	fullmsg = "| ** Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	if ropt:
		return "\n" + border + "\n" + fullmsg + "\n" + border + "\n";
	else:
		print("\n" + border + "\n" + fullmsg + "\n" + border + "\n");

#############################################################################

def getOutdir(indir, prefix, suffix, stime):
# Retrieves full input directory name and proper output directory name for other scripts.
	if not os.path.isdir(indir):
		errorOut(0, "-i must be a valid directory path");
		sys.exit();
	indir = os.path.abspath(indir);
	filelist = os.listdir(indir);
	used = [0];
	for each in filelist:
		if each.find("-" + prefix) != -1:
			used.append(int(each[:each.index("-")]));
	outdir = os.path.join(indir, str(max(used)+1) + "-" + prefix + "-" + stime + suffix);

	return indir, outdir;

#############################################################################

def spacedOut(string, totlen):
# Properly adds spaces to the end of a string to make it a given length
	spaces = " " * (totlen - len(string));
	return string + spaces;

#############################################################################

def filePrep(filename, header=""):
# Writes over a file, header optional (if no header just pass "")
	f = open(filename, "w");
	f.write(header);
	f.close();

#############################################################################
def listCheck(lst):
# Checks if all elements in a list are the same or not.
	if lst.count(lst[0]) == len(lst):
		return True;
	else:
		return False;

#############################################################################

def checkAlign(seqdict):
# Checks if a set of sequences is of equal length.
	aln_flag = True;
	seqlen = len(seqdict[list(seqdict.keys())[0]]);
	for title in seqdict:
		if len(seqdict[title]) != seqlen:
			aln_flag = False;
	return aln_flag;

#############################################################################

def defaultOutFile(input_name, file_flag, suffix="", output_init=False):
	i = 2;
	if suffix != "" and suffix[0] != "-":
		suffix = "-" + suffix;
	if not output_init:
		output, ext = list(os.path.splitext(input_name));
	# If the user did not specify an output file name, take the base of the input file name.
	else:
		output, ext = list(os.path.splitext(output_init));
	# Otherwise, use the user specified option.
	if not file_flag:
		if output[-1] in ["\\", "/"]:
			output = output[:-1];
		output = output + suffix + "-1" + ".txt";
	else:
		output = output + suffix + "-1" + ext;

	while os.path.exists(output) or os.path.exists(os.path.splitext(output)[0]):
		output = os.path.splitext(output);
		output = output[0][:output[0].rindex("-")+1] + str(i) + output[1];
		i += 1;
	# If the chosen output file exists, this will continually add 1 to a counter label at the end
	# of the file until a new file is chosen that does not exist.

	return output, (i-1);

#############################################################################

def defaultOutDir(input_name, file_flag, suffix="", output_init=False):
	i = 2;
	if suffix != "" and suffix[0] != "-":
		suffix = "-" + suffix;
	if not output_init:
		output = os.path.splitext(input_name)[0].rstrip("/").rstrip("\\") + suffix;
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

def getFileReader(i_name):
# Check if a file is gzipped, and if so set gzip as the file reader. Otherwise, read as a normal text file.
	try:
		gzip_check = gzip.open(i_name).read(1);
		reader = gzip.open;
	except:
		reader = open;
	return reader;

#############################################################################

def getFileLenRead(i_name):
# Gets the number of lines in a file.
	num_lines = 0;
	for line in getFileReader(i_name)(i_name): num_lines += 1;
	return float(num_lines);

#############################################################################

def dsum(*dicts):
# Given a list of dictionaries, this function merges them into a single dictionary, summing values of common keys.
    ret = defaultdict(int);
    for d in dicts:
        for k, v in d.items():
            ret[k] += v;
    return dict(ret);

#############################################################################

def chunks(l, n):
# Splits a list l into even chunks of size n.
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

#############################################################################

def runTime(msg=False):
	if msg:
		print("###### " + msg + " ######");
	print("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])))
	print("# Script call: " + " ".join(sys.argv))
	print("# Runtime: " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"));
	print("----------");

#############################################################################
##########################################################################################################################################################
#SEQUENCE FORMAT READERS AND WRITERS
##########################################################################################################################################################
#FASTA

def fastaGetLists(i_name):
#fastaGetLists reads a file and parses (separates) each FASTA sequence in the file into two corresponding lists:
#one containing the title line of the sequence and another containing the sequence itself. This function also
#returns the input and output file names.

	infile = open(i_name, "r");
	inseqs = infile.read();
	infile.close();
	#This block reads the input file.

	seqs = [''] * inseqs.count('>');
	titles = [''] * inseqs.count('>');

	k = 0;

	if len(seqs) > 1:
		for k in range(len(seqs) - 1):
			titles[k] = inseqs[:inseqs.index('\n')];

			seqs[k] = inseqs[inseqs.index('\n') + 1:inseqs.index('>', inseqs.index('\n')) - 1];
			tmp = inseqs.index('>', inseqs.index('\n'));
			inseqs = inseqs[tmp:];

			titles[k + 1] = inseqs[:inseqs.index('\n')];
			seqs[k + 1] = inseqs[inseqs.index('\n'):];

	else:
		titles[0] = inseqs[:inseqs.index('\n')];
		seqs[0] = inseqs[inseqs.index('\n') + 1:];

	for i in range(len(seqs)):
		seqs[i] = seqs[i].replace('\n', '');
	#The above lines parse the seqs file by creating a list with the correct amount of elements (one per gene) based on the count of FASTA title lines ('>'). It
	#then reads through the input string and stores all gene sequences in a separate list element. Finally, all newlines are removed.

	return titles, seqs;

#############################################################################

def fastaReader(i_name, meth="dict"):
# This function takes an input file and determines if it is a compressed or
# uncompressed FASTA file (.fa). If it doesn't end with .fa it returns it to
# be skipped.
	#try:
	if i_name.endswith(".fa.gz") or i_name.endswith(".fas.gz") or i_name.endswith(".faa.gz") or i_name.endswith(".fna.gz") or i_name.endswith(".fasta.gz"):
		if meth == "dict":
			seqs = fastaGetDictCompressed(i_name);
		elif meth == "ind":
			seqs = fastaReadInd(i_name);
	elif i_name.endswith(".fa") or i_name.endswith(".fas") or i_name.endswith(".faa") or i_name.endswith(".fna") or i_name.endswith(".fasta"):
		if meth == "dict":
			seqs = fastaGetDict(i_name);
		elif meth == "ind":
			seqs = fastaReadInd(i_name);
	else:
		return None, i_name;
	#except:
	#	return None, i_name;
	if seqs == {}:
		return None, i_name;
	return seqs, False;
	# If the reading of the file was successful, it returns the sequences. If not,
	# it returns the file name to be recorded as being skipped.

#############################################################################

def relabelHeader(title, new_label, header_delim, ropt):
	#print title, new_label, header_delim, ropt;
	if ropt == 1:
		new_title = ">" + new_label + header_delim + title[1:];
	elif ropt == 2:
		new_title = ">" + new_label;
	elif ropt == 3:
		new_title = title + header_delim + new_label;
	return new_title;

#############################################################################

def getOutFile(fasta_file, file_flag, out_dest, label):
	if file_flag:
		if out_dest:
			outfilename = out_dest;
		else:
			outfilename = os.path.splitext(fasta_file)[0] + "." + label + os.path.splitext(fasta_file)[1];
	else:
		outfilename = os.path.join(out_dest, os.path.splitext(os.path.basename(fasta_file))[0] + "." + label + os.path.splitext(fasta_file)[1]);
	return outfilename;

#############################################################################

def fastaGetDict(i_name):
#fastaGetDict reads a FASTA file and returns a dictionary containing all sequences in the file with 
#the key:value format as title:sequence.

	seqdict = {};
	for line in open(i_name, "r"):
		line = line.replace("\n", "");
		if line[:1] == '>':
			curkey = line;
			seqdict[curkey] = "";
		else:
			seqdict[curkey] = seqdict[curkey] + line;

	return seqdict;

#############################################################################

def fastaGetDictCompressed(i_name):
#fastaGetDict reads a FASTA file and returns a dictionary containing all sequences in the file with 
#the key:value format as title:sequence.
	import gzip

	seqdict = {};
	for line in gzip.open(i_name, "rb"):
		line = line.decode().replace("\n", "");
		if line[:1] == '>':
			curkey = line;
			seqdict[curkey] = "";
		else:
			seqdict[curkey] = seqdict[curkey] + line;

	return seqdict;


#############################################################################

def fastaGetFileInd(i_name):
#fastaGetFileInd reads a FASTA file and returns a dictionary containing file indexes for each title
#and sequence with the key:value format as [title start index]:[sequence start index]

	infile = open(i_name, "r");
	indList = [];
	firstflag = 0;
	curlist = [];

	line = "derp";

	while line != '':
		line = infile.readline();
		if line[:1] == '>':
			if firstflag == 1:
				curseqend = infile.tell() - len(line) - 1;
				curlist.append(curseqend);
				indList.append(curlist);
				curlist = [];

			curtitlestart = infile.tell() - len(line);
			curtitleend = infile.tell() - 1;
			curseqstart = infile.tell();

			curlist.append(curtitlestart);
			curlist.append(curtitleend);
			curlist.append(curseqstart);

			firstflag = 1;

	curseqend = infile.tell() - len(line) - 1;
	curlist.append(curseqend);
	indList.append(curlist);

	infile.close();
	return indList;
		
#############################################################################

def getFastafromInd(i_name, titlestart, titleend, seqstart, seqend):
#This takes the file index for a corresponding FASTA title and sequence (as retrieved by
#fastaGetFileInd and returns the actual text of the title and the sequence.

	infile = open(i_name, "r");

	infile.seek(titlestart);
	title = infile.read(titleend - titlestart);

	infile.seek(seqstart);
	seq = infile.read(seqend - seqstart);

	infile.close();

	title = title.replace("\n", "");
	seq = seq.replace("\n", "");

	return title, seq;

#############################################################################

def fastaReadInd(i_name):
# fastaGetFileInd reads a FASTA file and returns a dictionary containing file indexes for each title
# and sequence with the key:value format as [title start index]:[sequence start index]

	try:
		gzip_check = gzip.open(i_name).read(1);
		reader = gzip.open;
	except:
		reader = open;
	# Check if the genotype likelihood file is gzipped, and if so set gzip as the file reader. Otherwise, read as a normal text file.
	# MAKE SURE THIS WORKS WITH GZIPPED FILES!!

	with reader(i_name, "rb") as infile:
		fasta, first, curlist = {}, False, [];
		line = "derp";
		while line != '':
			line = infile.readline();
			if line[:1] == '>':
				if first:
					curseqend = infile.tell() - len(line) - 1;
					curlist.append(curseqend);
					fasta[cur_title] = curlist;
					curlist = [];

				#cur_title = line[1:].strip().split(" ")[0];
				cur_title = line.strip();
				curtitlestart = infile.tell() - len(line);
				curtitleend = infile.tell() - 1;
				curseqstart = infile.tell();

				curlist.append(curtitlestart);
				curlist.append(curtitleend);
				curlist.append(curseqstart);

				first = True;

		curseqend = infile.tell() - len(line);
		curlist.append(curseqend);
		fasta[cur_title] = curlist;

	return fasta;
		
#############################################################################

def fastaGetInd(i_name, inds):
# This takes the file index for a corresponding FASTA title and sequence (as retrieved by
# fastaGetFileInd and returns the actual text of the title and the sequence.

	titlestart, titleend, seqstart, seqend = inds;

	with open(i_name, "rb") as infile:
		infile.seek(titlestart);
		title = infile.read(titleend - titlestart);

		infile.seek(seqstart);
		seq = infile.read(seqend - seqstart);

	title = title.replace("\n", "");
	seq = seq.replace("\n", "");

	return title, seq;

#############################################################################

def writeFasta(seqdict, o_name):
#Writes a sequence dictionary in FASTA format to the given output file.
	for each in seqdict:
		title = each;
		seq = seqdict[each];
		seq = seq.replace(" ", "");
		outfile = open(o_name, "a");
		outfile.write(title);
		outfile.write("\n");
		a = 1;
		for base in seq:
			outfile.write(base);
			if a % 60 == 0:
				outfile.write("\n");
			a = a + 1;
		outfile.write("\n");
		outfile.close();

#############################################################################

def writeSeq(o_name, seq, title):
#A function to write the sequences to the output file in FASTA format. The
#sequences will be written 60 characters to a line.
	seq = seq.replace(" ", "");
	outfile = open(o_name, "a");
	outfile.write(title);
	outfile.write("\n");
	a = 1;
	for base in seq:
		outfile.write(base);
		if a % 60 == 0:
			outfile.write("\n");
		a = a + 1;
	outfile.write("\n");
	outfile.close();

#############################################################################

def writeSeqOL(o_name, seq, title):
#A funtion to write the sequences to the output file in FASTA format. The
#sequences will be written to one line.
	seq = seq.replace(" ", "");
	outfile = open(o_name, "a");
	outfile.write(title);
	outfile.write("\n");
	outfile.write(seq);
	outfile.write("\n");
	outfile.close();

##########################################################################################################################################################
#NEXUS

def nexusGetDict(i_name):
#A function to read sequences from a Nexus formatted file and store them in
#a dictionary.
	seqs = {};
	seqflag = 0;

	for line in open(i_name):
		if line.lower().find("begin data;") != -1:
			seqflag = 1;
			continue;

		if line.lower().find("end;") != -1:
			break;

		if seqflag == 1:
			if line[0] != "\t":
				tmpline = line.replace("\n","").split(" ");
				tmpline = [_f for _f in tmpline if _f];
				seqs[tmpline[0]] = tmpline[1];

	return seqs;

#############################################################################

def writeNexus(seqdict, o_name):
#A function which takes a sequence dictionary and output file name and writes
#the sequences to the file in Nexus format.
	outfile = open(o_name, "w");
	outfile.write("#NEXUS\n\nBegin data;\n");

	ntaxa = len(seqdict);
	nsites = len(seqdict[list(seqdict.keys())[0]]);

	outline = "\tDimensions ntax=" + str(ntaxa) + " nchar=" + str(nsites) + ";\n";
	outfile.write(outline);
	outline = "\tFormat datatype=protein gap=-;\n";
	outfile.write(outline);
	outline = "\tMatrix\n";
	outfile.write(outline);

	interval = len(max(list(seqdict.keys()), key=len)) + 3;

	for title in seqdict:
		newtitle = title.replace(" ",":");
		spaces = " " * (interval - len(newtitle));
		seq = seqdict[title];

		outline = newtitle + spaces + seq + "\n";
		outfile.write(outline);

	outline = "\t;\n";
	outfile.write(outline);
	outline = "End;";
	outfile.write(outline);

	outfile.close();

##########################################################################################################################################################
#PHYLIP

def phylipGetInterleaved(i_name):
#A function to read an interleaved (annoying) Phylip formatted file. Stores the
#sequences in a dictionary and also returns the info in the first line about
#seq len and num taxa.
	seqdict = {};

	i = 0;

	for line in open(i_name, "r"):
		line = line.replace("\n", "");
		if i == 0:
			firstline = line;
			i = i + 1;
			continue;

		if (line.count("A") + line.count("T") + line.count("C") + line.count("G") + line.count("-") + line.count("X") + line.count("N")) != len(line):
			curkey = line;
			seqdict[curkey] = "";
		else:
			seqdict[curkey] = seqdict[curkey] + line;

	return seqdict, firstline;

#############################################################################

def phylipGetDict(i_name):
#A function to read a Phylip formatted file. Returns the sequences in a
#dictionary format.
	seqdict = {};
	titlelist = [];

	i = 0;
	titlecount = 0;

	for line in open(i_name, "r"):
		if i == 0:
			firstline = line;
			i = i + 1;
			continue;
		if line == "" or line == "\n":
			continue;

		if line[:1] != " ":
			curtitle = line[:line.index(" ")];
			titlelist.append(curtitle);

			curseq = line[line.index(" "):];
			curseq = curseq.replace(" ","");
			curseq = curseq.replace("\n","");

			seqdict[curtitle] = curseq;
			continue;

		if titlecount >= len(titlelist):
			titlecount = 0;

		curseq = line.replace(" ","");
		curseq = curseq.replace("\n","");

		seqdict[titlelist[titlecount]] = seqdict[titlelist[titlecount]] + curseq;
		titlecount = titlecount + 1;
		
	return seqdict, firstline;

#############################################################################

def writePhylip(seqs, o_name):
#This function takes a sequence dictionary and writes it in Phylip format to
#the output file name.
	outfile = open(o_name, "w");
	outline = " " + str(len(seqs)) + " " + str(len(seqs[list(seqs.keys())[0]])) + "\n";
	outfile.write(outline);

	interval = len(max(list(seqs.keys()), key=len)) + 3;

	for title in seqs:
		newtitle = title.replace(" ", ":");
		spaces = " " * (interval - len(newtitle));
		seq = seqs[title];

		outline = newtitle + spaces + seq + "\n";
		outfile.write(outline);

	outfile.close();

#############################################################################
