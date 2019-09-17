
import math
from collections import defaultdict
def fqFunc(cur_reads):
    read_info = {
        'reads' : 0,
        'length' : 0,
        'num_sites' : 0, 
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
    
        read_info['num_sites'] += readlen;
        cur_bin = math.floor(readlen/5)*5;
        read_info['read_lens'][cur_bin] += 1;
        # Get the read length bins if that option is specified.
    
        for char in read['seq']:
            read_info['base_comp'][char] += 1;
        # Get the base counts if the base composition option is specified.
    
        for pos in range(readlen):
            read_info['qual_pos'][pos] += ord(read['qual'][pos]) - 33;
            read_info['site_pos'][pos] += 1;
        # Get the quality sums for each position and total number of sites for each position if the qual option is specified.
    
    return read_info;
    