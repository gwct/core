#############################################################################
# FASTQ functions -- used by fastq_stats.py
# Gregg Thomas
# September 2019
#############################################################################

import core, sys, os, math
from collections import defaultdict

#############################################################################

def printStats(num_reads, len_sum, genome_size, base_comp, read_lens, qual_pos, site_sum, pyv):

    if pyv == '2':
        bar = u'\u2588'.encode('utf-8');
    elif pyv == '3':
        bar = "█";

    print("TOTAL READS:         " + str(num_reads));
    avg_len = len_sum / num_reads;
    print("AVERAGE READ LENGTH: " + str(round(avg_len, 2)));
    if genome_size:
        coverage = (num_reads * avg_len) / genome_size;
        print("COVERAGE:            " + str(round(coverage, 2)));
    # Count the number of reads.

    #####
    print("\nBASE COMPOSITION:");
    base_sum = sum(base_comp.values());
    bases = ["A","T","C","G","N"];
    for base in bases:
        if base_comp[base] == 0:
            bc = "0.00";
        else:
            bc = str(round(base_comp[base] / base_sum, 2))

        base_comp[base] = bc;

    print("".join([b + "     " for b in bases]));
    outline = "";
    for base in bases:
        outline = outline + base_comp[base] + "  ";
    print(outline);
    # Calculate base composition
    #####

    #####
    print("\nRough read length distribution:");
    print("Read length....Proportion of reads");

    b = min(read_lens) - 5;
    max_bin = max(read_lens);

    print("                       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9       1.0")
    print("               ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|")
    while b <= max_bin+5:
        b_str = str(b);
        while len(b_str) != 3:
            b_str = "0" + b_str;
        if b in read_lens:
            read_lens[b] = round(read_lens[b] / num_reads, 2) * 100;
            if read_lens[b] == 0.00:
                print(b_str + "............|"); 
            else:
                print(b_str + "............" + bar*int(read_lens[b]));
        else:
            print(b_str + "............");
        b += 5;
    print("               ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|")
    print("                       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9       1.0")
    # Read length proportion histogram.
    #####

    #####
    # Average quality per position plot
    print("\nAverage quality score per read position:");
    print("Position....Quality");
    print("           0        10        20        30        40");
    print("           |---------|---------|---------|---------|-");
    for pos in sorted(qual_pos.keys()):
        qual_pos[pos] = round(qual_pos[pos] / site_sum[pos]);

        pos_str = str(pos+1);
        while len(pos_str) < 3:
            pos_str = "0" + pos_str;

        if qual_pos[pos] == 0:
            print(pos_str + "........!");
        else:
            print(pos_str + "........|" + "-" * (qual_pos[pos]-1) + chr(qual_pos[pos] + 33));
    print("           |---------|---------|---------|---------|-");
    print("           0        10        20        30        40");
    # Average quality per position plot
    #####

####################

def countReads(fastq_files, genome_size, dirflag, paired):
# This function counts the number of sequences and positions in a given set of FASTA files.
    py_vers = sys.version[0];

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
            print("\n--------------------\n")
            print("FILE:\t" + fqf)

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

            rotator, divisor = 0,1000;
            for line in core.getFileReader(fqf)(fqf):
                line = line.decode().strip();
                line_counter += 1;
                if line.startswith("@"):
                    line_counter = 0;

                if line_counter == 1:
                    num_reads += 1;
                    len_sum += len(line);
                    
                    cur_bin = math.floor(len(line)/5)*5;
                    read_lens[cur_bin] += 1;
                    
                    if paired and fastq_file.index(fqf) == 0:
                        total_reads += 1;
                        total_len_sum += len(line);
                        total_read_lens[cur_bin] += 1;

                    for char in line:
                        base_comp[char] += 1;
                        if paired and fastq_file.index(fqf) == 0:
                            total_base_comp[char] += 1;

                if line_counter == 3:
                    for pos in range(len(line)):
                        qual = ord(line[pos]) - 33;

                        qual_pos[pos] += qual;
                        total_qual_pos[pos] += qual;
                        site_sum[pos] += 1;
                        total_site_sum[pos] += 1;

            printStats(num_reads, len_sum, genome_size, base_comp, read_lens, qual_pos, site_sum, py_vers);

    if dirflag:
        print("\n=======================================================================");
        print("\nCOMBINED STATS FOR ALL FILES:")
        printStats(total_reads, total_len_sum, genome_size, total_base_comp, total_read_lens, total_qual_pos, total_site_sum, py_vers);



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