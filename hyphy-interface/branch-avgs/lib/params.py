#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#
# This dictionary will also be used to track output.
#############################################################################

import sys
import timeit
import lib.core as CORE

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
    globs = {
        'version' : 'Beta 1.0',
        'releasedate' : "October 2021",
        'authors' : "Emily Kopania, Gregg Thomas",
        'doi' : '',
        'http' : '',
        'github' : '',
        'starttime' : timeit.default_timer(),
        'startdatetime' : CORE.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

        'tree-info-file' : False,
        'csv-rate-dir' : False,
        'filter-file' : False,
        'subset-file' : False,
        'outdir' : False,
        'output-file' : False,
        # Input options

        'csv-files' : [],
        'filter-files' : [],
        'subset-files' : [],
        # Gene lists

        'branches' : {},
        # The main dictionary for the tree info

        'num-branches' : 0,

        'clade-index' : 1,
        # The index of the 'clade' column in the tree info csv
        # TODO: This could just be read from the headers...

        'headers' : False,
        # The headers in the input file + new

        'run-name' : 'branch-avgs',
        'logfilename' : 'branch-avgs.errlog',
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

    globs['logfilename'] = "branch-avgs-" + globs['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    #globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################