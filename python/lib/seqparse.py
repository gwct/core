#############################################################################
# Functions to read and parse sequences.
# Will split others out of core over time.
# 04.2020
#############################################################################

import sys
import core

#############################################################################

def dnaCheck(seq, i_name):
# dnaCheck does several things to ensure the input sequence is in proper format (DNA FASTA).

    dna_symbols = ["A","T","C","G","N","-"];

    seq = seq.upper();
    # It first converts the sequence to all uppercase letters.
    seq = seq.replace('U', 'T');
    # The function then converts any U's to T's, changing RNA to DNA.

    for n in range(len(seq)):
    # This loop then reads through every letter of the sequence and makes sure they are all nucleotide symbols.
    # If they are not, an error is printed and the program exits.
        if seq[n] not in dna_symbols:        
            print("\nError! Input sequence from the file", i_name, "is not a nucleotide sequence. Please input only DNA sequences.\n")
            sys.exit();
    return seq;
    # If the input sequence checks out, it is returned to the call.

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

def premStopCheck(seq, frame=1, allowlastcodon=False, rmlast=False):
# Checks a coding sequence for premature stop codons. Default is first frame.
    stop_codons = ["TAG", "TAA", "TGA", "UAG", "UAA", "UGA"];
    seq = seq.upper();

    if frame not in [1,2,3]:
        sys.exit(" * SEQ ERROR: premStopCheck: Invalid reading frame input: " + str(frame));

    if frame == 2:
        seq = seq[1:];
    if frame == 3:
        seq = seq[2:];

    codon_list = [ seq[i:i+3] for i in range(0, len(seq), 3) ];
    #codon_list_orig = codon_list.copy();
    codon_list_orig = [ codon for codon in codon_list ];
    while codon_list[-1] == "---":
        codon_list = codon_list[:-1]

    is_stop = False;
    for c in range(len(codon_list)):
        if codon_list[c] in stop_codons:
            if c+1 == len(codon_list):
                if rmlast:
                    codon_list_orig[c] = "NNN"; 
                if not allowlastcodon:
                    is_stop = True;
            else:
                is_stop = True;

    return is_stop, "".join(codon_list_orig);

#############################################################################

def ntToCodon(nt_seq):
# Splits a nucleotide sequence into a list of codons
    #assert len(nt_seq) % 3 == 0, "\nOUT OF FRAME NUCLEOTIDE SEQUENCE! " + str(len(nt_seq));
    codon_seq = [(nt_seq[i:i+3]) for i in range(0, len(nt_seq), 3)];
    return codon_seq;

#############################################################################

def revComp(seq, iupac=True):
# Returns the reverse complement of a nucleotide sequence.
    if iupac:
        complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y' : 'R', 'R' : 'Y', 'S' : 'S', 'W' : 'W', 'K' : 'M', 'M' : 'K', 'B' : 'V', 'V' : 'B', 'D' : 'H', 'H' : 'D', 'N' : 'N' };
    else:
        complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'C', 'N' : 'N' };

    return "".join(complement.get(base, base) for base in reversed(seq.upper()));

#############################################################################

def bioTranslator(i_name):
# This function takes a FILE of DNA sequences in FASTA format and translates it
# to the corresponding AA sequence. It then returns that sequence as a dictionary.
#
# This function may become deprecated as the call to fastaGetDict can be performed
# before the call to this function, increasing functionality. See newbioTranslator
# below.

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

#############################################################################

def newbioTranslator(seq):
# This function takes a DNA sequence as a single string and returns the
# corresponding AA sequence.

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

def yabt(codon_seq):
# Yet another biotranslator with a more logical way to store the genetic code

    standard_code = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
            'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W',
            '---':'-', 'NNN':'N'
        }
    # From: https://www.geeksforgeeks.org/dna-protein-python-3/

    aa_seq = "";
    for codon in codon_seq:
        codon = codon.upper();
        if len(codon) != 3:
            aa_seq += "X";
        elif "N" in codon:
            aa_seq += "X";
        else:
            aa_seq += standard_code[codon];

    return aa_seq;

    #aa_seq = [ "N" if "N" in codon else standard_code[codon] for codon in codon_seq ];

    #return "".join(aa_seq);

#############################################################################

##########################################################################################################################################################
#SEQUENCE FORMAT READERS AND WRITERS
##########################################################################################################################################################
#FASTA

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
    #    return None, i_name;
    if seqs == {}:
        return None, i_name;
    return seqs, False;
    # If the reading of the file was successful, it returns the sequences. If not,
    # it returns the file name to be recorded as being skipped.

#############################################################################

def fastaGetLists(i_name):
# fastaGetLists reads a file and parses (separates) each FASTA sequence in the file into two corresponding lists:
# one containing the title line of the sequence and another containing the sequence itself. This function also
# returns the input and output file names.

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
    # The above lines parse the seqs file by creating a list with the correct amount of elements (one per gene) based on the count of FASTA title lines ('>'). It
    # then reads through the input string and stores all gene sequences in a separate list element. Finally, all newlines are removed.

    return titles, seqs;

#############################################################################

def fastaGetDict(i_name):
# fastaGetDict reads a FASTA file and returns a dictionary containing all sequences in the file with 
# the key:value format as title:sequence.

    seqdict = {};
    for line in open(i_name, "r"):
        if line == "\n":
            continue;
        line = line.replace("\n", "");
        if line[0] == '>':
            curkey = line;
            seqdict[curkey] = "";
        else:
            seqdict[curkey] = seqdict[curkey] + line;

    return seqdict;

#############################################################################

def fastaGetDictCompressed(i_name):
# fastaGetDict reads a FASTA file and returns a dictionary containing all sequences in the file with 
# the key:value format as title:sequence.
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

def fastaReadSeqs(filename, header_sep=False):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 
# Returns dictionary with the key:value format as title:sequence.

    import gzip
    from itertools import groupby

    compression = core.detectCompression(filename);

    if compression == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif compression == "none":
        file_stream = open(filename); 
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
        readstr = lambda s : s.strip();
    # Read the lines of the file depending on the compression level
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level -- for compressed files we also need to decode
    # each string in the iterators below.

    seqdict = {};
    # A dictionary of sequences:
    # <sequence id/header> : <sequence>

    for header_obj in fa_iter:
        #print(header_obj)
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        curkey = header[1:];
        # This removes the ">" character from the header string to act as the key in seqdict

        if header_sep:
            curkey = curkey.split(header_sep)[0];

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        seqdict[curkey] = seq;
        # Save the sequence in seqdict

    return seqdict;

#############################################################################

def fastaReadInd(filename):
# readFastaInd reads a FASTA file and returns a list of lists for each entry in the file
# Each sublist is [start index of title, end index of title, start index of seq, end index of seq]

    #import gzip
    # For possible future support

    index_list = [];
    # The list of indices to return

    compression = core.detectCompression(filename);
    # Detects the compression of the file
    # For now, just a check to make sure it is uncompressed

    if compression == "gz":
        sys.exit("COMPRESSED FILES CURRENTLY NOT SUPPORTED IN fastaReadInd().");
        #infile = gzip.open(filename);        
    elif compression == "none":
        infile = open(filename); 
    # Open the file or error out based on compression level

    cur_list = [];
    first = True;
    # Initialize the first list and the flag
    # The flag is so we don't try to add a sequence the first time we encounter
    # a header

    line = "init";
    # Initialize the line string so the loop starts

    while line != '':
        line = infile.readline();
        # Read each line in the file, one by one, until the end of the file

        if line[:1] == ">":
        # We add things to the cur_list when we encounter the header char

            if not first:
            # Every time we encounter a header char except the first, we
            # need to update the main list and re-initialize the cur_list

                cur_seq_end = infile.tell() - len(line) - 1;
                cur_list.append(cur_seq_end);
                # Get the end of the sequence as the end of the previous line

                index_list.append(cur_list);
                cur_list = [];
                # Add the current list to the main and reset for the next sequence

            cur_title_start = infile.tell() - len(line);
            # Get the position of the beginning of the current line

            cur_title_end = infile.tell() - 1;
            # Get the position of the end of the current lint

            cur_seq_start = infile.tell();
            # Get the position of the start of th sequence

            cur_list.append(cur_title_start);
            cur_list.append(cur_title_end);
            cur_list.append(cur_seq_start);
            # Add all of these to the cur_list

            first = False;
            # Be sure to set the flag to False after the first sequence 

    cur_seq_end = infile.tell() - len(line) - 1;
    cur_list.append(cur_seq_end);
    index_list.append(cur_list);
    # At the end, we need to add in the last sequence

    infile.close();
    # Close the file

    return index_list;
        
#############################################################################

def getSeqfromInd(filename, titlestart, titleend, seqstart, seqend, compression="none"):
# This takes the file index for a corresponding FASTA title and sequence (as retrieved by
# fastaReadInd and returns the actual text of the title and the sequence.

    #import gzip;
    # For possible future support

    if compression == "gz":
        sys.exit("COMPRESSED FILES CURRENTLY NOT SUPPORTED IN getSeqfromInd().");
        #infile = gzip.open(filename);
    elif compression == "none":
        infile = open(filename);
    # Open the file or error out based on compression level

    infile.seek(titlestart);
    title = infile.read(titleend - titlestart);
    # Get the title string by moving to the start and reading for the length of the title

    infile.seek(seqstart);
    seq = infile.read(seqend - seqstart);
    # Get the sequence string by moving to the start and reading the length of the sequence

    infile.close();
    # Close the file

    title = title.strip();
    seq = seq.strip();
    # Get rid of newlines

    return title, seq;

#############################################################################

def writeSeq(o_name, seq, title):
# A function to write the sequences to the output file in FASTA format. The
# sequences will be written 60 characters to a line.
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
# A funtion to write the sequences to the output file in FASTA format. The
# sequences will be written to one line.
    seq = seq.replace(" ", "");
    outfile = open(o_name, "a");
    outfile.write(title);
    outfile.write("\n");
    outfile.write(seq);
    outfile.write("\n");
    outfile.close();

#############################################################################

##########################################################################################################################################################
#NEXUS - not up to date

def nexusGetDict(i_name):
# A function to read sequences from a Nexus formatted file and store them in
# a dictionary.
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
# A function which takes a sequence dictionary and output file name and writes
# the sequences to the file in Nexus format.
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
#PHYLIP - not up to date

def phylipGetInterleaved(i_name):
# A function to read an interleaved (annoying) Phylip formatted file. Stores the
# sequences in a dictionary and also returns the info in the first line about
# seq len and num taxa.
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
# A function to read a Phylip formatted file. Returns the sequences in a
# dictionary format.
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
# This function takes a sequence dictionary and writes it in Phylip format to
# the output file name.
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
