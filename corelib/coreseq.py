############################################################
# Functions to read and parse sequences.
# Will split others out of core over time.
# 04.2020
############################################################

import sys

############################################################

def premStopCheck(seq, frame=1, allowlastcodon=False, rmlast=False):
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
                if allowlastcodon:
                    is_stop = False;
                else:
                    is_stop = True;
            else:
                is_stop = True;

    return is_stop, "".join(codon_list_orig);

############################################################