#############################################################################
# FASTQ functions -- used by fastq_stats.py
# Gregg Thomas
# September 2019
#############################################################################

import core, sys, os, math, gzip, copy
from collections import defaultdict
import multiprocessing as mp
import fqfunc as fq

#############################################################################

def textOut(stats, globs):
# Outputs fastq stats in a nice format to read.

    if globs['pyv'] == '2':
        bar = u'\u2588'.encode('utf-8');
    elif globs['pyv'] == '3':
        bar = u'\u2588';
    pad = 25;

    core.PW("\n" + core.spacedOut("FILE:", pad) + stats['file'], globs['outtxt']);

    if globs['reads']:
        core.PW(core.spacedOut("TOTAL READS:", pad) + str(stats['reads']), globs['outtxt']);
    if globs['lens']:
        avg_len = stats['sites'] / stats['reads'];
        core.PW(core.spacedOut("AVERAGE READ LENGTH:", pad) + str(round(avg_len, 2)), globs['outtxt']);
    if globs['reads'] and globs['lens'] and globs['genome_size']:
        coverage = (stats['reads'] * avg_len) / globs['genome_size'];
        core.PW(core.spacedOut("COVERAGE:", pad) + str(round(coverage, 2)), globs['outtxt']);
    # Count the number of reads.

    #####
    if globs['bc']:
        core.PW("\nBASE COMPOSITION:", globs['outtxt']);
        base_sum = sum(stats['base_comp'].values());
        bases = ["A","T","C","G","N"];
        for base in bases:
            if stats['base_comp'][base] == 0:
                bc = "0.00";
            else:
                bc = str(round(stats['base_comp'][base] / base_sum, 3))
                while len(bc) < 4:
                    bc += "0";

            stats['base_comp'][base] = bc;

        core.PW("".join([b + "      " for b in bases]), globs['outtxt']);
        outline = "";
        for base in bases:
            outline = outline + stats['base_comp'][base] + "  ";
        core.PW(outline, globs['outtxt']);
    # Calculate base composition
    #####

    #####
    if globs['lens']:
        core.PW("\nRough read length distribution:", globs['outtxt']);
        core.PW("Read length....Proportion of reads", globs['outtxt']);

        b = min(stats['read_lens']) - 5;
        max_bin = max(stats['read_lens']);

        core.PW("                       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9       1.0", globs['outtxt'])
        core.PW("               ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|", globs['outtxt'])
        while b <= max_bin+5:
            b_str = str(b);
            while len(b_str) != 3:
                b_str = "0" + b_str;
            if b in stats['read_lens']:
                stats['read_lens'][b] = round(stats['read_lens'][b] / stats['reads'], 2) * 100;
                if stats['read_lens'][b] == 0.00:
                    core.PW(b_str + "............|", globs['outtxt']); 
                else:
                    core.PW(b_str + "............" + bar*int(stats['read_lens'][b]), globs['outtxt']);
            else:
                core.PW(b_str + "............", globs['outtxt']);
            b += 5;
        core.PW("               ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|", globs['outtxt'])
        core.PW("                       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9       1.0", globs['outtxt'])
        core.PW("Read length....Proportion of reads", globs['outtxt']);
    # Read length proportion histogram.
    #####

    #####
    # Average quality per position plot
    if globs['qual']:
        core.PW("\nAverage quality score per read position:", globs['outtxt']);
        core.PW("Position...Quality", globs['outtxt']);
        core.PW("           0        10        20        30        40", globs['outtxt']);
        core.PW("           |---------|---------|---------|---------|-", globs['outtxt']);
        for pos in sorted(stats['qual_pos'].keys()):
            stats['qual_pos'][pos] = round(stats['qual_pos'][pos] / stats['site_pos'][pos]);

            pos_str = str(pos+1);
            while len(pos_str) < 3:
                pos_str = "0" + pos_str;

            if stats['qual_pos'][pos] == 0:
                core.PW(pos_str + "........!", globs['outtxt']);
            else:
                core.PW(pos_str + "........|" + "-" * (stats['qual_pos'][pos]-1) + chr(stats['qual_pos'][pos] + 33), globs['outtxt']);
        core.PW("           |---------|---------|---------|---------|-", globs['outtxt']);
        core.PW("           0        10        20        30        40", globs['outtxt']);
        core.PW("Position...Quality", globs['outtxt']);
    # Average quality per position plot
    #####

####################

def csvOut(stats, globs):
# Outputs fastq stats in a nice format to parse.

    with open(globs['outcsv'], 'a') as csvfile:
        outline = [stats['file']];

        if globs['reads']:
            outline.append(str(stats['reads']));
        # Num reads output

        if globs['lens']:
            avg_len = stats['sites'] / stats['reads'];
            outline.append(str(round(avg_len, 3)));
        # Read length output

        if globs['reads'] and globs['lens'] and globs['genome_size']:
            coverage = round((stats['reads'] * avg_len) / globs['genome_size'], 3);
            outline.append(str(coverage));
        # Genome coverage output

        if globs['bc']:
            base_sum = sum(stats['base_comp'].values());
            bases = ["A","T","C","G","N"];
            for base in bases:
                if stats['base_comp'][base] == 0:
                    bc = "0.00";
                else:
                    bc = str(round(stats['base_comp'][base] / base_sum, 3))
                    while len(bc) < 4:
                        bc += "0";

                stats['base_comp'][base] = str(bc);

            outline.extend([stats['base_comp'][base] for base in bases]);
        # Base composition output

        if globs['qual']:
            qual_sum = sum(stats['qual_pos'].values());
            avg_qual = round(qual_sum / stats['sites'], 3);
            # Compute the average quality of all bases.

            sorted_qual_pos = sorted(stats['qual_pos'].keys())
            qual_pos_first = { pos : stats['qual_pos'][pos] for pos in sorted_qual_pos[:5] }
            qual_pos_last = { pos : stats['qual_pos'][pos] for pos in sorted_qual_pos[-5:] }

            sorted_count_pos = sorted(stats['site_pos'].keys())
            count_pos_first = { pos : stats['site_pos'][pos] for pos in sorted_count_pos[:5] }
            count_pos_last = { pos : stats['site_pos'][pos] for pos in sorted_count_pos[-5:] }            

            qual_sum_first = sum(qual_pos_first.values());
            count_sum_first = sum(count_pos_first.values());
            avg_qual_first = round(qual_sum_first / count_sum_first, 3);
            # Comput the average quality of the first 5 bases in all reads.

            qual_sum_last = sum(qual_pos_last.values());
            count_sum_last = sum(count_pos_last.values());
            avg_qual_last = qual_sum_last / count_sum_last;
            # Compute the average quality of the last 5 bases of the longest reads. This doesn't work well because reads are of different lengths.

            outline.extend([str(avg_qual), str(avg_qual_first)]);
        # Quality output

        csvfile.write(",".join(outline) + "\n");
        

####################

def fqProcessRead(lines=None):
# Takes a list of four lines assumed to be the four lines from a fastq read and stores them in a dictionary.
    ks = ['title', 'seq', 'optional', 'qual'];
    return {k: v for k, v in zip(ks, lines)};

####################

def chunks(l, n):
# Splits a list l into even chunks of size n.
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

####################

def dsum(*dicts):
    # Given a list of dictionaries, this function merges them into a single dictionary, summing values of common keys.
    ret = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            ret[k] += v
    return dict(ret)

####################
def combineResults(ocounts, fcounts, r):
    ocounts['reads'] += r['reads'];
    ocounts['sites'] += r['sites'];
    ocounts['read_lens'] = dsum(ocounts['read_lens'], r['read_lens']);
    ocounts['base_comp'] = dsum(ocounts['base_comp'], r['base_comp']);
    ocounts['qual_pos'] = dsum(ocounts['qual_pos'], r['qual_pos']);
    ocounts['site_pos'] = dsum(ocounts['site_pos'], r['site_pos']);

    fcounts['reads'] += r['reads'];
    fcounts['sites'] += r['sites'];
    fcounts['read_lens'] = dsum(fcounts['read_lens'], r['read_lens']);
    fcounts['base_comp'] = dsum(fcounts['base_comp'], r['base_comp']);
    fcounts['qual_pos'] = dsum(fcounts['qual_pos'], r['qual_pos']);
    fcounts['site_pos'] = dsum(fcounts['site_pos'], r['site_pos']);

    return ocounts, fcounts

####################
def createFunc(globs):
# This function dynamically writes a function depending on the input options. This prevents having an if statement
# for each option for each line in a fastq file. Pretty cool...

    func = '''
import math
from collections import defaultdict
def fqFunc(cur_reads):
    read_info = {
        'reads' : 0,
        'length' : 0,
        'sites' : 0, 
        'read_lens' : defaultdict(int),
        'base_comp' : { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 },
        'qual_pos' : defaultdict(int),
        'site_pos' : defaultdict(int)
        }
    # All the info possibily collected for a set of reads is compiled in this dictionary.

    for read in cur_reads:
        read_info['reads'] += 1;
        readlen = len(read['seq']);
        # Always get the number of reads and the read lengths.
    '''

    if globs['lens']:
        func += '''
        read_info['sites'] += readlen;
        cur_bin = math.floor(readlen/5)*5;
        read_info['read_lens'][cur_bin] += 1;
        # Get the read length bins if that option is specified.
    '''

    if globs['bc']:
        func += '''
        for char in read['seq']:
            read_info['base_comp'][char] += 1;
        # Get the base counts if the base composition option is specified.
    '''

    if globs['qual']:
        func += '''
        for pos in range(readlen):
            read_info['qual_pos'][pos] += ord(read['qual'][pos]) - 33;
            read_info['site_pos'][pos] += 1;
        # Get the quality sums for each position and total number of sites for each position if the qual option is specified.
    '''

    func += '''
    return read_info;
    '''

    funcpath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "fqfunc.py");
    with open(funcpath, "w") as funcfile:
        funcfile.write(func);
    # Write the function to a file. This file has been imported above.

####################

def countReads(fastq_files, globs):
# This function counts the number of sequences and positions in a given set of FASTA files.

    createFunc(globs);

    overall_counts = {
        'file' : globs['dirflag'],
        'reads' : 0,
        'sites' : 0,
        'read_lens' : defaultdict(int),
        'base_comp' : { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 },
        'qual_pos' : defaultdict(int),
        'site_pos' : defaultdict(int)
    }
    # total_reads: The total number of reads parsed in all input files
    # total_sites: The total number of sites in all input files
    # total_read_lens: The count of read lengths for every read in every file
    # total_base_comp: The base compositions for every file
    # total_qual_pos: The sum of quality scores from all positions in all reads in all files
    # total_site_pos: The total number of sites for all positions in all reads in all files

    for fastq_file in fastq_files:
        for fqf in fastq_file:
            print("\n# --------------------\n");
            print("# " + core.getTime() + " READING FILE :\t" + fqf);

            file_counts = {
                'file' : fqf,
                'reads' : 0,
                'sites' : 0,
                'read_lens' : defaultdict(int),
                'base_comp' : { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 },
                'qual_pos' : defaultdict(int),
                'site_pos' : defaultdict(int)
            }
            # num_reads: The total number of reads parse in the current file
            # len_sum: The sum length of all reads in the current file
            # site_sum: The total number of sites at each position in all reads
            # qual_pos: The sum of quality scores for every position in every read
            # read_lens: The number of reads of every length in the current file
            # base_comp: The base compositions for the current file

            if fqf.endswith(".gz"):
                lread = lambda l : l.decode().rstrip();
            else:
                lread = lambda l : l.rstrip();
            # Declare a lambda function to parse lines in the file based on the file type. This prevents having an if statement for each line read.

            lines, reads = [], [];
            reads_per_chunk = 10000;
            pool = mp.Pool(processes=globs['procs']);

            print("# " + core.getTime() + " Processing reads in chunks of size " + str(reads_per_chunk) + " using " + str(globs['procs']) + " processes.");
            i, i_start = 0, 1;
            for line in core.getFileReader(fqf)(fqf):
                lines.append(lread(line));
                if len(lines) == 4:
                    i += 1;
                    read = fqProcessRead(lines);
                    reads.append(read);
                    lines = [];

                if len(reads) == globs['procs']*reads_per_chunk:
                    print("# " + core.getTime() + " Processing reads " + str(i_start) + "-" + str(i));
                    i_start = i + 1;
                    read_chunks = list(chunks(reads, reads_per_chunk));
                    reads = [];
                    for result in pool.imap(fq.fqFunc, (read_chunk for read_chunk in read_chunks)):
                        overall_counts, file_counts = combineResults(overall_counts, file_counts, result);

            if reads != []:
                print("# " + core.getTime() + " Processing reads " + str(i_start) + "-" + str(i));
                read_chunks = list(chunks(reads, reads_per_chunk));
                reads = [];
                for result in pool.imap(fq.fqFunc, (read_chunk for read_chunk in read_chunks)):
                    overall_counts, file_counts = combineResults(overall_counts, file_counts, result);

            print("# " + core.getTime() + " Done parsing...");
            if not globs['summary']:
                if globs['outtxt']:
                    print("# " + core.getTime() + " Printing stats and writing to text file...");
                    textOut(copy.deepcopy(file_counts), globs);
                if globs['outcsv']:
                    print("\n# " + core.getTime() + " Writing stats to csv file...");
                    csvOut(file_counts, globs);
            # Output the stats for the current file.
 
    print("\n# --------------------\n");
    if globs['outtxt']:
        print("# " + core.getTime() + " Printing TOTAL stats and writing to text file...");
        textOut(copy.deepcopy(overall_counts), globs);
    if globs['outcsv']:
        print("\n# " + core.getTime() + " Writing TOTAL stats to csv file...");
        csvOut(overall_counts, globs);
    # Output summary stats for ALL files combined.

    print("# " + core.getTime() + " Done!");
    print("\n# =======================================================================");
