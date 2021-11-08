#!/usr/bin/env python3
#############################################################################
# A template for python scripts to make things easier to start
#############################################################################

import sys
import os
import multiprocessing as mp
import lib.core as CORE
import lib.params as params
import lib.opt_parse as OP
# import lib.other as OTHER

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# branch_avgs version " + globs['version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.

    print("#");
    print("# " + "=" * 125);
    #print(CORE.welcome());
    #if "-h" not in sys.argv:
    #    print("            A description of the program\n");
    # A welcome banner.

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ##########################

    ## DO STUFF HERE

    step = "A step";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    things_done = 0;

    step_start_time = CORE.report_step(globs, step, step_start_time, "SUCCESS: " + str(len(things_done)) + " things done!");
    # Status update

    ## DO STUFF HERE

    ##########################

    ## DO SOMETHING PARALLEL

    step = "Parallel step";
    step_start_time = CORE.report_step(globs, step, False, "Processed 0 / " + str(globs['total']) + " things...", full_update=True);
    # Status update

    with mp.Pool(processes=globs['num-procs']) as pool:
        counter = 1;
        new_data = {};
        for result in pool.imap_unordered(BRANCHES.branchSum, ((globs['branches'][branch], branch, globs['csv-rate-dir'], globs['filter-files'], globs['subset-files']) for branch in branches)):
            core.PWS("# " + core.getDateTime() + " Finished branch " + str(counter) + " / " + branch_num);
            #print(result[0], result[1]);
            
            new_branches[result[0]] = result[1];
            #branches.update(result);
            cur_branch_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(branch_num) + " / " + str(globs['total']) + " things...", full_update=True);
            branch_num += 1;

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=True);
    # Status update

    ## DO SOMETHING PARALLEL

    ##########################

    CORE.endProg(globs);

#############################################################################

