#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#
# This dictionary will also be used to track output.
#############################################################################

import sys
import timeit
import lib.core as PC

#############################################################################

class StrictDict(dict):
# This prevents additional keys from being added to the global params dict in
# any other part of the code, just to help me limit it's scope
# https://stackoverflow.com/questions/32258706/how-to-prevent-key-creation-through-dkey-val
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError("{} is not a legal key of this StricDict".format(repr(key)));
        dict.__setitem__(self, key, value);

#############################################################################

def init():
    globs_init = {
        'version' : 'Beta 1.0',
        'releasedate' : "October 2021",
        'authors' : "Gregg Thomas",
        'doi' : '',
        'http' : '',
        'github' : '',
        'starttime' : timeit.default_timer(),
        'startdatetime' : PC.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

        'gxf-file' : False,
        'fa-file' : False,
        'gxf-type' : False,
        # Input with annotation file and genome file

        'longest-isoform-ids' : [],
        # The list of longest isoform IDs

        'in-seq' : False,
        'in-seq-type' : False,
        # Input by a directory with many fasta files or a single multi-fasta

        'gxf-compression' : 'none',
        'seq-compression' : 'none',
        # The type of compression used for input sequence files

        'annotation' : {},
        # The main annotation storage dict. A nested structure for genes, transcripts, and coding exons.
        # <gene id> : { <header>, <start coord>, <end coord>, <strand>, 
        #                   { <transcript id> : <transcript header>, <transcript start coord>, <transcript end coord> <transcript strand>, 
        #                       { <exon id> : <exon header>, <exon start coord>, <exon end coord>, <exon strand> } } }

        'seqs' : {},
        # The dictionary in which to read the input sequences

        'out-header-prefix' : "",
        'splitter' : False,
        'id-repl' : False,

        'header' : False,
        'out-header' : 'retain',
        # I'm thinking this could be "retain" or "id", along with an option to add a species id
        'get-prot-ids' : False,
        'tid-to-pid' : {},
        # User options

        'outdir' : '',
        'out-id-file' : '',
        'out-seq-file' : '', 
        'tid-to-pid-file' : '',
        'run-name' : 'isofilter',
        'logfilename' : 'isofilter.errlog',
        'logdir' : '',
        'overwrite' : False,
        # I/O options

        'num-procs' : 1,
        # Number of processes to use

        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'qstats' : False,
        'norun' : False,
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "isofilter-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################