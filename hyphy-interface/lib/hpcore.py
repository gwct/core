############################################################
# Support functions for the PAML parser.
# 11.2020
############################################################

import sys, datetime

############################################################

def readMeta(metafile):
    features, first = {}, True;
    for line in open(metafile):
        if line[0] == "#":
            continue;
        if first:
            first = False;
            continue;
        line = line.strip().split("\t");
        fid, ttype, chrome, start, end, strand, gid, num_cds = line;
        features[fid] = { 'gid' : gid, 'chrome' : chrome, 'start' : start, 'end' : end, 'strand' : strand, 'num-cds' : int(num_cds) };
    return features;

############################################################

def runTime(msg=False, writeout=False):
	if msg:
		if not msg.startswith("#"):
			msg = "# " + msg;
		PWS(msg, writeout);

	PWS("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])), writeout)
	PWS("# Script call:    " + " ".join(sys.argv), writeout)
	PWS("# Runtime:        " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"), writeout);
	PWS("# ----------------", writeout);

############################################################

def PWS(o_line, o_stream=False, std_stream=True):
# Function to print a string AND write it to the file.
	if std_stream:
		print(o_line);
	if o_stream:
		o_stream.write(o_line + "\n");

############################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a string to make it a given length
	spaces = sep * (totlen - len(string));
	return string + spaces;

############################################################

def getDateTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

############################################################

        # if scaff == "ScmyWZ3_7747_HRSCAF_7900":
        #     if int(end) < 35225145:
        #         scaff = "ScmyWZ3_7747_HRSCAF_7900_R";
        #         chrome = "chrX_R";
        #     elif int(start) >= 35225145:
        #         scaff = "ScmyWZ3_7747_HRSCAF_7900_NR"
        #         chrome = "chrX_NR";
        #     else:
        #         chrome = "chrX";
        # elif scaff in scaff_to_chr:
        #     chrome =  scaff_to_chr[scaff];
        # else:
        #     chrome = "NA";