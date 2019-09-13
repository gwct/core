#############################################################################
# FASTQ functions -- used by fastq_stats.py
# Gregg Thomas
# September 2019
#############################################################################

import core, sys, os, math
from collections import defaultdict
from Bio import SeqIO as seq

#############################################################################

def printStats(stats, globs):

    if globs['pyv'] == '2':
        bar = u'\u2588'.encode('utf-8');
    elif globs['pyv'] == '3':
        bar = u'\u2588';

    if globs['reads']:
        print("TOTAL READS:         " + str(stats['num_reads']));
    if globs['lens']:
        avg_len = stats['len_sum'] / stats['num_reads'];
        print("AVERAGE READ LENGTH: " + str(round(avg_len, 2)));
    if globs['reads'] and globs['lens'] and globs['genome_size']:
        coverage = (stats['num_reads'] * avg_len) / globs['genome_size'];
        print("COVERAGE:            " + str(round(coverage, 2)));
    # Count the number of reads.

    #####
    if globs['bc']:
        print("\nBASE COMPOSITION:");
        base_sum = sum(stats['base_comp'].values());
        bases = ["A","T","C","G","N"];
        for base in bases:
            if stats['base_comp'][base] == 0:
                bc = "0.00";
            else:
                bc = str(round(stats['base_comp'][base] / base_sum, 2))
                while len(bc) != 4:
                    bc += "0";

            stats['base_comp'][base] = bc;

        print("".join([b + "     " for b in bases]));
        outline = "";
        for base in bases:
            outline = outline + stats['base_comp'][base] + "  ";
        print(outline);
    # Calculate base composition
    #####

    #####
    if globs['lens']:
        print("\nRough read length distribution:");
        print("Read length....Proportion of reads");

        b = min(stats['read_lens']) - 5;
        max_bin = max(stats['read_lens']);

        print("                       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9       1.0")
        print("               ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|")
        while b <= max_bin+5:
            b_str = str(b);
            while len(b_str) != 3:
                b_str = "0" + b_str;
            if b in stats['read_lens']:
                stats['read_lens'][b] = round(stats['read_lens'][b] / stats['num_reads'], 2) * 100;
                if stats['read_lens'][b] == 0.00:
                    print(b_str + "............|"); 
                else:
                    print(b_str + "............" + bar*int(stats['read_lens'][b]));
            else:
                print(b_str + "............");
            b += 5;
        print("               ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|")
        print("                       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9       1.0")
    # Read length proportion histogram.
    #####

    #####
    # Average quality per position plot
    if globs['qual']:
        print("\nAverage quality score per read position:");
        print("Position....Quality");
        print("           0        10        20        30        40");
        print("           |---------|---------|---------|---------|-");
        for pos in sorted(stats['qual_pos'].keys()):
            stats['qual_pos'][pos] = round(stats['qual_pos'][pos] / stats['site_sum'][pos]);

            pos_str = str(pos+1);
            while len(pos_str) < 3:
                pos_str = "0" + pos_str;

            if stats['qual_pos'][pos] == 0:
                print(pos_str + "........!");
            else:
                print(pos_str + "........|" + "-" * (stats['qual_pos'][pos]-1) + chr(stats['qual_pos'][pos] + 33));
        print("           |---------|---------|---------|---------|-");
        print("           0        10        20        30        40");
    # Average quality per position plot
    #####

####################

def fqProcessRead(lines=None):
    ks = ['title', 'seq', 'optional', 'qual'];
    return {k: v for k, v in zip(ks, lines)};

####################

def countReads(fastq_files, globs):
# This function counts the number of sequences and positions in a given set of FASTA files.

    total_reads, total_len_sum, total_sites = 0,0,0;

    total_site_sum = defaultdict(int);
    total_qual_pos = defaultdict(int);
    
    total_read_lens = defaultdict(int);
    total_base_comp = { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 };
    # total_reads: The total number of reads parsed in all input files
    # total_len_sum: The sum of the lengths of all reads from all input files
    # total_sites: The total number of sites in all input files
    # total_site_sum: The total number of sites for all positions in all reads in all files
    # total_qual_pos: The sum of quality scores from all positions in all reads in all files
    # total_read_lens: The count of read lengths for every read in every file
    # total_base_comp: The base compositions for every file

    for fastq_file in fastq_files:
        for fqf in fastq_file:
            print("\n--------------------\n");
            print("FILE :\t" + fqf);
            #numlines = core.getFileLen(fqf);
            #print("LINES:\t" + str(numlines));

            num_reads, len_sum, line_counter, num_sites = 0,0,0,0;

            site_sum = defaultdict(int);
            qual_pos = defaultdict(int);

            read_lens = defaultdict(int);
            base_comp = { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 };
            # num_reads: The total number of reads parse in the current file
            # len_sum: The sum length of all reads in the current file
            # line_counter: To tell which part of the fastq entry we're on
            # num_sites: The total number of sites parsed in the current file
            # site_sum: The total number of sites at each position in all reads
            # qual_pos: The sum of quality scores for every position in every read
            # read_lens: The number of reads of every length in the current file
            # base_comp: The base compositions for the current file

            fq_lines = core.getFileReader(fqf)(fqf).read().decode().split("\n");
            lines = [];
            reads = [];
            for line in (fq_lines):
                lines.append(line.rstrip());
                if len(lines) == 4:
                    read = fqProcessRead(lines);
                    reads.append(read);
                    lines = [];

            print("READ READS");
            print(len(reads));
            sys.exit();                
            
            
            #fqhandler = core.getFileReader(fqf)(fqf, "r");

            numbars, donepercent, i = 0,[],0;
            #with core.getFileReader(fqf)(fqf) as infile:
            lines = [];
            for line in seq.parse(fqhandler, "fastq"):
                numbars, donepercent = core.loadingBar(i, numlines, donepercent, numbars);
                i += 1;
                print(line);
                sys.exit();


                #lines.append(line);
                lines.append(line.decode().rstrip());
                if len(lines) == 4:
                    read = fqProcessRead(lines);

                    if globs['reads']:
                        num_reads += 1;
                    
                    if globs['lens']:
                        readlen = len(read['seq'])
                        len_sum += readlen;                           
                        cur_bin = math.floor(readlen/5)*5;
                        read_lens[cur_bin] += 1;
                        
                    if globs['paired'] and fastq_file.index(fqf) == 0:
                        if globs['reads']:
                            total_reads += 1;
                        if globs['lens']:
                            total_len_sum += readlen;
                            total_read_lens[cur_bin] += 1;

                    if globs['bc']:
                        for char in read['seq']:
                            base_comp[char] += 1;
                            if globs['paired'] and fastq_file.index(fqf) == 0:
                                total_base_comp[char] += 1;

                    if globs['qual']:
                        for pos in range(readlen):
                            qual = ord(read['qual'][pos]) - 33;
                            qual_pos[pos] += qual;
                            total_qual_pos[pos] += qual;
                            site_sum[pos] += 1;
                            total_site_sum[pos] += 1;

                    lines = []

            pstring = "100.0% complete.";
            sys.stderr.write('\b' * len(pstring) + pstring + "\n");

            stats = {
                'num_reads' : num_reads,
                'len_sum' : len_sum,
                'base_comp' : base_comp,
                'read_lens' : read_lens,
                'qual_pos' : qual_pos,
                'site_sum' : site_sum
            }

            printStats(stats, globs);


    if globs['dirflag']:
        stats = {
            'num_reads' : total_reads,
            'len_sum' : total_len_sum,
            'base_comp' : total_base_comp,
            'read_lens' : total_read_lens,
            'qual_pos' : total_qual_pos,
            'site_sum' : total_site_sum
        }

        print("\n=======================================================================");
        print("\nCOMBINED STATS FOR ALL FILES:")
        printStats(stats, globs);



	# 	seqs, skip = core.fastaReader(fasta_file);
	# 	if skip:
	# 		fa_skip.append(fasta_file);
	# 		continue;
	# 	total_files += 1;
	# 	for title in seqs:
	# 		total_seq += 1;
	# 		if disp_file == 1:
	# 			print(title + "\t" + len(seqs[title]));
	# 		total_pos += len(seqs[title]);

	# print("\n" + core.getTime() + " Done!");
	# print("-----");
	# print("Total FASTA files:\t", total_files);
	# print("Total sequences:\t", total_seq);
	# print("Total positions:\t", total_pos);
	# if fa_skip != []:
	# 	print("The following", str(len(fa_skip)), "file(s) were skipped because they couldn't be read as fasta files: ", ",".join([os.path.basename(f) for f in fa_skip]));
	# print("=======================================================================");

#############################################################################