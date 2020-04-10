############################################################
# Functions to read and parse sequences.
# Will split others out of core over time.
# 04.2020
############################################################

import sys

############################################################

def premStopCheck(seq, frame=1):
    stop_codons = ["TAG", "TAA", "TGA", "UAG", "UAA", "UGA"];
    seq = seq.upper();

    if frame not in [1,2,3]:
        sys.exit(" * SEQ ERROR: premStopCheck: Invalid reading frame input: " + str(frame));

    if frame == 2:
        seq = seq[1:];
    if frame == 3:
        seq = seq[2:];

    codon_list = [ seq[i:i+3] for i in range(0, len(seq), 3) ];
    for c in range(len(codon_list)):
        if c+1 == len(codon_list):
            continue;

        if codon_list[c] in stop_codons:
            return True;

    return False

############################################################