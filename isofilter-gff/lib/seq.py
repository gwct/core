#############################################################################
# Functions to read sequences for degenotate
#############################################################################

import sys
import os
import gzip
import lib.core as CORE
from itertools import groupby

############################################################################# 

def readFasta(filename, globs):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 

    if globs['seq-compression'] == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif globs['seq-compression'] == "none":
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

        #curkey = header[1:header.index(" ")];
        curkey = header[1:];
        # This removes the ">" character from the header string to act as the key in seqdict
        # TODO: Need to decide if splitting the input sequence header on " " is good for most cases, or if we should
        #       add this as a user option?

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        seqdict[curkey] = seq;
        # Save the sequence in the dictionary

    file_stream.close();

    return seqdict;

#############################################################################

def readSeqs(globs):
# A function that reads an entire genome fasta file into memory

    step = "Detecting compression of FASTA file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['seq-compression'] = CORE.detectCompression(globs['fa-file']);
    if globs['seq-compression'] == "none":
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No compression detected");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + globs['seq-compression'] + " detected");
    # Detect the compression of the input sequence file

    step = "Reading FASTA file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['seqs'] = readFasta(globs['fa-file'], globs);

    if len(globs['seqs']) == 0:
        step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(len(globs['seqs'])) + " seqs read");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['seqs'])) + " seqs read");
    # Read the input sequence file

    #print(list(globs['genome-seqs'].keys()))

    return globs;

#############################################################################

