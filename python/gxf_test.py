import core
import seqparse as seq
import sys
import gxfparse as gxf

# infilename = "../test-data/Mus_musculus.GRCm38.pep.all.fa"
# infilename = "../test-data/Mus_musculus.GRCm38.cds.all.fa.gz"


# seqs = seq.fastaReadInd(infilename);
# print(len(seqs));

# for s in seqs:
#     print(s);
#     t, seq = seq.getSeqfromInd(infilename, s[0], s[1], s[2], s[3]);
#     print(t);
#     print(seq);
#     sys.exit();


infilename = "../test-data/Mus_musculus.GRCm38.99.gff3.gz";

annotation, compression, tid_to_gid = gxf.getExons(infilename, coding_only=True);