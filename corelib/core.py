#############################################################################
#CORE functions
#Gregg Thomas
#August 2013-present
#############################################################################

import string
import sys
import re
import subprocess
import datetime

#############################################################################

def dnaCheck(seq, i_name):
#dnaCheck does several things to ensure the input sequence is in proper format (DNA FASTA).

	dna_symbols = ["A","T","C","G","N","-"];

	seq = string.upper(seq);
	#It first converts the sequence to all uppercase letters.
	seq = seq.replace('U', 'T');
	#The function then converts any U's to T's, changing RNA to DNA.

	for n in xrange(len(seq)):
	#This loop then reads through every letter of the sequence and makes sure they are all nucleotide symbols.
	#If they are not, an error is printed and the program exits.
		if seq[n] not in dna_symbols:		
			print "\nError! Input sequence from the file", i_name, "is not a nucleotide sequence. Please input only DNA sequences.\n"
			sys.exit();
	return seq;
	#If the input sequence checks out, it is returned to the call.

#############################################################################

def loadingBar(counter, length, done, bars):
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

	percent = float(counter) / float(length) * 100.0;
	percentdone = int(percent);

	p = str(percent)
	pstring = " " + p[:5] + "% complete.";

	if percentdone % 2 == 0 and done != None and percentdone not in done:
		loading = "";
		loading = "[";
		j = 0;
		while j <= bars:
			loading = loading + "*";
			j = j + 1;
		while j < 50:
			loading = loading + "-";
			j = j + 1;
		loading = loading + "]";

		loading = loading + "                 ";
		sys.stderr.write('\b' * len(loading) + loading);

		done.append(percentdone);
		bars = bars + 1;

	sys.stderr.write('\b' * len(pstring) + pstring);

	return bars, done;

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
	return datetime.datetime.now().strftime("%m.%d.%Y-%I:%M:%S");

#############################################################################

def printWrite(o_name, o_line):
#Function to print a string AND write it to the file.
	print o_line;
	f = open(o_name, "a");
	f.write(o_line + "\n");
	f.close();

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
				tmpline = filter(None, tmpline);
				seqs[tmpline[0]] = tmpline[1];

	return seqs;

#############################################################################

def writeNexus(seqdict, o_name):
#A function which takes a sequence dictionary and output file name and writes
#the sequences to the file in Nexus format.
	outfile = open(o_name, "w");
	outfile.write("#NEXUS\n\nBegin data;\n");

	ntaxa = len(seqdict);
	nsites = len(seqdict[seqdict.keys()[0]]);

	outline = "\tDimensions ntax=" + str(ntaxa) + " nchar=" + str(nsites) + ";\n";
	outfile.write(outline);
	outline = "\tFormat datatype=protein gap=-;\n";
	outfile.write(outline);
	outline = "\tMatrix\n";
	outfile.write(outline);

	interval = len(max(seqdict.keys(), key=len)) + 3;

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
	outline = " " + str(len(seqs)) + " " + str(len(seqs[seqs.keys()[0]])) + "\n";
	outfile.write(outline);

	interval = len(max(seqs.keys(), key=len)) + 3;

	for title in seqs:
		newtitle = title.replace(" ", ":");
		spaces = " " * (interval - len(newtitle));
		seq = seqs[title];

		outline = newtitle + spaces + seq + "\n";
		outfile.write(outline);

	outfile.close();

#############################################################################
