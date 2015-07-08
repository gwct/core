#############################################################################
#Core functions
#Gregg Thomas
#August 2013
#############################################################################


import string
import sys
import re
import subprocess
import datetime

#############################################################################

def dnaCheck(seq, iName):
#dnaCheck does several things to ensure the input sequence is in proper format (DNA FASTA).

	seq = string.upper(seq);
	#It first converts the sequence to all uppercase letters.

	seq = seq.replace('U', 'T');
	#The function then converts any U's to T's, changing RNA to DNA.

	for n in range(len(seq)):
	#This loop then reads through every letter of the sequence and makes sure they are all nucleotide symbols.
	#If they are not, an error is printed and the program exits.

		if seq[n] != 'A' and seq[n] != 'T' and seq[n] != 'C' and seq[n] != 'G' and seq[n] != '-':
			
			print "\nError! Input sequence from the file", iName, "is not a nucleotide sequence. Please input only DNA sequences.\n"
			sys.exit();
	return seq;
	#If the input sequence checks out, it is returned to the call.

#############################################################################

def fastaGetLists(inputFileName):
#fastaGetLists reads a file and parses (separates) each FASTA sequence in the file into two corresponding lists:
#one containing the title line of the sequence and another containing the sequence itself. This function also
#returns the input and output file names.

	inFile = open(inputFileName, "r");
	inputSeqs = inFile.read();
	inFile.close();
	#This block reads the input file.

	seqs = [''] * inputSeqs.count('>');
	titleLines = [''] * inputSeqs.count('>');

	k = 0;

	if len(seqs) > 1:
		for k in range(len(seqs) - 1):

			titleLines[k] = inputSeqs[:inputSeqs.index('\n')];

			seqs[k] = inputSeqs[inputSeqs.index('\n') + 1:inputSeqs.index('>', inputSeqs.index('\n')) - 1];
			tmp = inputSeqs.index('>', inputSeqs.index('\n'));
			inputSeqs = inputSeqs[tmp:];


			titleLines[k + 1] = inputSeqs[:inputSeqs.index('\n')];
			seqs[k + 1] = inputSeqs[inputSeqs.index('\n'):];

	else:
		titleLines[0] = inputSeqs[:inputSeqs.index('\n')];
		seqs[0] = inputSeqs[inputSeqs.index('\n') + 1:];

	for i in range(len(seqs)):
		seqs[i] = seqs[i].replace('\n', '');
	#The above lines parse the seqs file by creating a list with the correct amount of elements (one per gene) based on the count of FASTA title lines ('>'). It
	#then reads through the input string and stores all gene sequences in a separate list element. Finally, all newlines are removed.


	#for q in range(len(seqs)):
	#	seqs[q] = dnaCheck(seqs[q], inputFileName);
	#Each sequence is passed to dnaCheck which parses the sequence and ensures it contains only nucleotide symbols.

	return titleLines, seqs;

#############################################################################

def fastaGetDict(inputFileName):
#fastaGetDicts reads a FASTA file and returns a dictionary containing all sequences in the file with 
#the key:value format as title:sequence.

	seqDict = {};

	for line in open(inputFileName, "r"):

		line = line.replace("\n", "");

		if line[:1] == '>':
			curkey = line;
			seqDict[curkey] = "";

		else:
			seqDict[curkey] = seqDict[curkey] + line;

	return seqDict;

#############################################################################

def fastaGetFileInd(inputFileName):
#fastaGetFileInd reads a FASTA file and returns a dictionary containing file indexes for each title
#and sequence with the key:value format as [title start index]:[sequence start index]

	inFile = open(inputFileName, "r");
	indList = [];
	firstflag = 0;
	curlist = [];

	line = "derp";

	while line != '':

		line = inFile.readline();

		if line[:1] == '>':

			if firstflag == 1:
				curseqend = inFile.tell() - len(line) - 1;
				curlist.append(curseqend);
				indList.append(curlist);
				curlist = [];

			curtitlestart = inFile.tell() - len(line);
			curtitleend = inFile.tell() - 1;
			curseqstart = inFile.tell();

			curlist.append(curtitlestart);
			curlist.append(curtitleend);
			curlist.append(curseqstart);

			firstflag = 1;

	curseqend = inFile.tell() - len(line) - 1;
	curlist.append(curseqend);
	indList.append(curlist);

	inFile.close();
	return indList;
		
#############################################################################

def getFastafromInd(inputFileName, titlestart, titleend, seqstart, seqend):
#This takes the file index for a corresponding FASTA title and sequence (as retrieved by
#fastaGetFileInd and returns the actual text of the title and the sequence.

	inFile = open(inputFileName, "r");

	inFile.seek(titlestart);
	title = inFile.read(titleend - titlestart);

	inFile.seek(seqstart);
	seq = inFile.read(seqend - seqstart);

	inFile.close();

	title = title.replace("\n", "");
	seq = seq.replace("\n", "");

	return title, seq;

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

def bioTranslator(inFileName):

	codons = [''] * 5;
	codons[0] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	codons[1] = "---M---------------M---------------M----------------------------";
	codons[2] = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	codons[3] = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	codons[4] = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

	dnaSeqs = fastaGetDict(inFileName);
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

def phylipGetDict(inputFileName):
#fastaGetDicts reads a FASTA file and returns a dictionary containing all sequences in the file with 
#the key:value format as title:sequence.

	seqDict = {};

	i = 0;

	for line in open(inputFileName, "r"):

		line = line.replace("\n", "");

		if i == 0:
			firstline = line;
			i = i + 1;
			continue;

		if (line.count("A") + line.count("T") + line.count("C") + line.count("G") + line.count("-") + line.count("X") + line.count("N")) != len(line):
			curkey = line;
			seqDict[curkey] = "";

		else:
			seqDict[curkey] = seqDict[curkey] + line;

	return seqDict, firstline;

#############################################################################

def phylipGetDict2(inputFileName):
#fastaGetDicts reads a FASTA file and returns a dictionary containing all sequences in the file with 
#the key:value format as title:sequence.

	seqDict = {};
	titleList = [];

	i = 0;
	titlecount = 0;

	for line in open(inputFileName, "r"):

		#line = line.replace("\n", "");

		if i == 0:
			firstline = line;
			i = i + 1;
			continue;
		if line == "" or line == "\n":
			continue;

		if line[:1] != " ":
			curtitle = line[:line.index(" ")];
			titleList.append(curtitle);

			curseq = line[line.index(" "):];
			curseq = curseq.replace(" ","");
			curseq = curseq.replace("\n","");

			seqDict[curtitle] = curseq;
			continue;

		if titlecount >= len(titleList):
			titlecount = 0;


		curseq = line.replace(" ","");
		curseq = curseq.replace("\n","");

		seqDict[titleList[titlecount]] = seqDict[titleList[titlecount]] + curseq;
		titlecount = titlecount + 1;
		

	return seqDict, firstline;

#############################################################################

def writeSeq(outputFileName, seq, title):
#A function to write the sequences to the output file in FASTA format. The
#sequences will be written 60 characters to a line.

	seq = seq.replace(" ", "");
	outFile = open(outputFileName, "a");
	outFile.write(title);
	outFile.write("\n");
	a = 1;
	for base in seq:
		outFile.write(base);
		if a % 60 == 0:
			outFile.write("\n");
		a = a + 1;
	outFile.write("\n");
	outFile.close();

#############################################################################

def writeSeqOL(outputFileName, seq, title):
#A funtion to write the sequences to the output file in FASTA format. The
#sequences will be written to one line.

	seq = seq.replace(" ", "");
	outFile = open(outputFileName, "a");
	outFile.write(title);
	outFile.write("\n");
	outFile.write(seq);
	outFile.write("\n");
	outFile.close();

#############################################################################

def variance(data):

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

def getFileLen(filename):
	p = subprocess.Popen(['wc', '-l', filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	result, err = p.communicate();
	if p.returncode != 0:
		raise IOError(err);
	return int(result.strip().split()[0]);

#############################################################################

def getDateTime():
	return datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

#############################################################################

def getTime():
	return datetime.datetime.now().strftime("%I:%M:%S");

#############################################################################

def getLogTime():
	return datetime.datetime.now().strftime("%m.%d.%Y-%I:%M:%S");

#############################################################################

def printWrite(filename, oline):
	print oline;
	f = open(filename, "a");
	f.write(oline + "\n");
	f.close();

#############################################################################
