#############################################################################
# Output functions for isofilter.
#############################################################################

import sys
import lib.core as CORE

#############################################################################

def writeLongestIsoformID(globs):
    step = "Writing IDs for longest isoforms";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");

    num_written = 0;
    num_prot_coding = 0;
    with open(globs['out-id-file'], "w") as outfile:
        for gene in globs['annotation']:
            if not globs['annotation'][gene]['transcripts']:
                continue;
            # Check if the current gene has at least 1 transcript, skip if not

            sorted_transcripts = sorted(globs['annotation'][gene]['transcripts'].items(), key = lambda k: k[1]['length'], reverse=True);
            # Sort the transcripts by length using this function to sort nested dictionaries
            # https://www.geeksforgeeks.org/python-sort-nested-dictionary-by-key/

            longest_transcript = sorted_transcripts[0][0];
            # The sorted function returns a list of tuples. Get the first value in the first tuple, which is the 
            # transcript ID with the longest length

            # print(sorted_transcripts);
            # for i in sorted_transcripts:
            #     print(i[0], i[1]['length']);
            # print(longest_transcript);
            # For debugging

            outfile.write(longest_transcript + "\n");
            num_written += 1;
            # Write the ID of the longest transcript

            if globs['get-prot-ids'] and longest_transcript in globs['tid-to-pid']:
                longest_transcript = globs['tid-to-pid'][longest_transcript];
                num_prot_coding += 1;
            # If the --prot option is enabled, convert between transcript ID and protein ID for protein coding transcripts here
            elif globs['id-repl'] and longest_transcript in globs['tid-to-pid']:
                longest_transcript = globs['tid-to-pid'][longest_transcript];
            # If the -repl option is enabled, convert between the original transcript ID and the edited one here
            

            globs['longest-isoform-ids'].append(longest_transcript);
            # Save the ID for the longest isoform to the global list

            #print("-----------\n");   

    if num_written == 0:
        step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(num_written) + " IDs written");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_written) + " IDs written");
    # Status update

    if globs['get-prot-ids']:
        if num_prot_coding == 0:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: " + str(num_prot_coding) + " protein coding transcripts found");
        else:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(num_prot_coding) + " protein coding transcripts found");
    # Info update if --prot is enabled 

#############################################################################

def writeLongestIsoformSeq(globs):
    step = "Writing sequences for longest isoforms";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");

    num_written = 0;
    with open(globs['out-seq-file'], "w") as outfile:
        for header in globs['seqs']:
            #print(header);
            for isoform_id in globs['longest-isoform-ids']:
                if isoform_id in header:
                    outfile.write(">" + globs['out-header-prefix'] + header + "\n");
                    outfile.write(globs['seqs'][header] + "\n");
                    num_written += 1;

    if num_written == 0:
        step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(num_written) + " sequences written");
    else:    
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_written) + " sequences written");

#############################################################################