
import math
from collections import defaultdict
def fqFunc(read):
    read_info = {
        'length' : 0, 
        'read_lens' : defaultdict(int),
        'base_comp' : { "A" : 0, "T" : 0, "C" : 0, "G" : 0, "N" : 0 },
        'qual_pos' : defaultdict(int),
        'site_pos' : defaultdict(int)
        }
    readlen = len(read['seq']);
    
    read_info['length'] = readlen;
    cur_bin = math.floor(readlen/5)*5;
    read_info['read_lens'][cur_bin] += 1;
    
    for char in read['seq']:
        read_info['base_comp'][char] += 1;
    
    for pos in range(readlen):
        read_info['qual_pos'][pos] += ord(read['qual'][pos]) - 33;
        read_info['site_pos'][pos] += 1;
    
    return read_info;
    