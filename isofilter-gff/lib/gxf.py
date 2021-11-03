#############################################################################
# Functions to handle gtf and gff annotation files.
#############################################################################

import sys
import os
import gzip
import lib.core as CORE

#############################################################################

def checkIDs(l, info, id_list, step, globs):
# Each time we read feature info and look for IDs, this checks to make sure only one
# ID is found. Probably unnecessary, but easy enough to check.

    if len(id_list) != 1:
        print("\n\n");
        print(id_list);
        print(info);
        print(l);
        print("\n\n");
        CORE.errorOut("GXF1", "Invalid number of IDs found during " + step, globs);

#############################################################################

# def readFeatures(feature_list, get_id=True, id_format, get_parent=True, parent_id_format, info_field_splitter)

#     return globs, num_read, num_missing
## TODO: Generalize the loops for each feature so we don't have all this copied code

#############################################################################

def read(globs):

    step = "Detecting compression of annotation file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['gxf-compression'] = CORE.detectCompression(globs['gxf-file']);
    if globs['gxf-compression'] == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No compression detected");
    else:
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + globs['gxf-compression'] + " detected");
    # Detect the compression of the input annotation file

    if globs['gxf-type'] == "gff":
        field_splitter = ";";
        gene_id_format = "ID=";
        transcript_id_format = "ID=";
        exon_id_format = "ID=";
        transcript_parent_format = "Parent=";
        exon_parent_format = "Parent=";

    elif globs['gxf-type'] == "gtf":
        field_splitter = "; ";
        gene_id_format = "gene_id ";
        transcript_id_format = "transcript_id ";
        exon_id_format = "exon_id ";
        transcript_parent_format = "gene_id ";
        exon_parent_format = "transcript_id ";
    # These outline the differences between GFF and GTF

    ####################

    step = "Reading genes";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");

    for line in reader(globs['gxf-file']):
        line = readline(line);
        # Read and parse the current line of the GXF file

        if "##FASTA" in line[0]:
            break;
        # Maker GFF files sometimes include the sequence of all transcripts at the end. We need to stop reading the file
        # at that point

        if line[0][0] == "#":
            continue;
        # Header/comment lines should be skipped. Note that this must come after the check for "##FASTA" above, or else
        # the file will keep being read into the sequences and error out

        feature_type, seq_header, start, end, strand, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[6], line[8].split(field_splitter);
        # Unpack the pertinent information from the current line into more readable variables

        if feature_type == "gene":
            feature_id = [ info_field for info_field in feature_info if info_field.startswith(gene_id_format) ];
            # Get the feature ID as a list of fields with the "ID=" prefix

            checkIDs(line, feature_info, feature_id, "gene id parsing", globs);
            # A quick check to make sure we have read only one ID

            feature_id = feature_id[0].replace(gene_id_format, "").replace("\"", "");
            # Unpack and parse the ID

            globs['annotation'][feature_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand, 'transcripts' : {} };
            # Add the ID and related info to the annotation dict. This includes an empty dict for transcripts to be stored in a similar way

    if len(globs['annotation']) == 0:
        step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(len(globs['annotation'])) + " genes read");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['annotation'])) + " genes read");
    # Status update

    # Read genes
    ####################

    step = "Reading transcripts";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");

    num_transcripts = 0;
    # A count of the number of transcripts read

    tid_to_gid = {};
    # A dictionary to link each transcript to it's gene ID, necessary for reading in exons in the next step
    # <transcript id> : <gene id>

    missing_genes = [];
    # Sometimes, transcripts will be present in the GXF file that have genes that don't have their own entry. This list
    # keeps track of those genes to report later

    for line in reader(globs['gxf-file']):
        line = readline(line);
        # Read and parse the current line of the GXF file

        if "##FASTA" in line[0]:
            break;
        # Maker GFF files sometimes include the sequence of all transcripts at the end. We need to stop reading the file
        # at that point.

        if line[0][0] == "#":
            continue;
        # Header/comment lines should be skipped. Note that this must come after the check for "##FASTA" above, or else
        # the file will keep being read into the sequences and error out.

        feature_type, seq_header, start, end, strand, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[6], line[8].split(field_splitter);
        # Unpack the pertinent information from the current line into more readable variables.

        if feature_type in ["transcript", "mRNA"]:
        # Skipping any 'unconfirmed_transcript'
            feature_id = [ info_field for info_field in feature_info if info_field.startswith(transcript_id_format) ];
            # Get the feature ID as a list of fields with the "ID=" prefix

            # print(feature_info);
            # print(feature_id);
            # sys.exit();

            checkIDs(line, feature_info, feature_id, "transcript id parsing", globs);
            # A quick check to make sure we have read only one ID

            feature_id = feature_id[0].replace(transcript_id_format, "").replace("\"", "");
            # Unpack and parse the ID

            parent_id = [ info_field for info_field in feature_info if info_field.startswith(transcript_parent_format) ];
            # Get the gene ID associated with the transcript as a list of fields with the "Parent=" prefix

            checkIDs(line, feature_info, parent_id, "transcript parent id parsing", globs);
            # A quick check to make sure we have read only one ID

            parent_id = parent_id[0].replace(transcript_parent_format, "").replace("\"", "");
            # Unpack and parse the gene ID

            if parent_id not in globs['annotation']:
                if parent_id not in missing_genes:
                    missing_genes.append(parent_id);
                continue;
            # If the gene doesn't have it's own entry in the GXF file, add it to the list of missing genes and
            # skip this transcript
            else:
                tid_to_gid[feature_id] = parent_id;
            # Otherwise, add the transcript id/gene id pair to the lookup dict

            transcript_length = end - start;

            # if globs['id-repl']:
            #     feature_id = feature_id.replace(globs['id-repl'][0], globs['id-repl'][1]);
            # If an ID replacement is specified in the options, do that replacement here SWITCHED TO BELOW

            globs['annotation'][parent_id]['transcripts'][feature_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'length' : transcript_length, 'strand' : strand, 'protein-id' : "" };
            # Add the ID and related info to the annotation dict. This includes an empty dict for exons to be stored in a similar way

            num_transcripts += 1;
            # Add to the number of transcripts read

    if num_transcripts == 0:
        step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(num_transcripts) + " transcripts read");
    else: 
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_transcripts) + " transcripts read");
    # Status update

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(len(missing_genes)) + " missing genes (present as parent of other feature, but not as its own feature)");  
    #print("# INFO: " + str(len(missing_genes)) + " missing genes (present as parent of other feature, but not as its own feature)");
    # Report the number of missing genes

    # Read transcripts
    ####################

    if globs['get-prot-ids']:

        step = "Getting protein IDs from coding exons";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");

        num_prot_ids = 0;
        # A count of the number of transcripts read

        missing_transcripts = [];
        # Sometimes, exons will be present in the GXF file that have transcripts that don't have their own entry. This list
        # keeps track of those transcripts to report later

        done = [];

        for line in reader(globs['gxf-file']):
            line = readline(line);
            # Read and parse the current line of the GXF file

            if "##FASTA" in line[0]:
                break;
            # Maker GFF files sometimes include the sequence of all transcripts at the end. We need to stop reading the file
            # at that point.

            if line[0][0] == "#":
                continue;
            # Header/comment lines should be skipped. Note that this must come after the check for "##FASTA" above, or else
            # the file will keep being read into the sequences and error out.

            feature_type, seq_header, start, end, strand, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[6], line[8].split(field_splitter);
            # Unpack the pertinent information from the current line into more readable variables.

            if feature_type in ["CDS"]:
                feature_id = [ info_field for info_field in feature_info if info_field.startswith(exon_id_format) ];
                # Get the feature ID as a list of fields with the "ID=" prefix

                #checkIDs(line, feature_info, feature_id, "exon id parsing", globs);
                # A quick check to make sure we have read only one ID

                if feature_id:
                    feature_id = feature_id[0].replace(exon_id_format, "").replace("\"", "");
                # Unpack and parse the ID

                parent_id = [ info_field for info_field in feature_info if info_field.startswith(exon_parent_format) ];
                # Get the transcript ID associated with the exon as a list of fields with the "Parent=" prefix

                checkIDs(line, feature_info, parent_id, "exon parent id parsing", globs);
                # A quick check to make sure we have read only one ID

                parent_id = parent_id[0].replace(exon_parent_format, "").replace("\"", "");
                # Unpack and parse the gene ID

                if parent_id not in tid_to_gid:
                    if parent_id not in missing_transcripts:
                        missing_transcripts.append(parent_id);
                    continue;
                # If the transcript doesn't have it's own entry in the GXF file, add it to the list of missing transcripts and
                # skip this exon
                elif parent_id in done:
                    continue;
                else:
                    gene_id = tid_to_gid[parent_id];
                # Otherwise, lookup the gene ID associated with the transcript to use below

                done.append(parent_id);

                protein_id = [ info_field for info_field in feature_info if info_field.startswith("protein_id=") ];
                # Get the transcript ID associated with the exon as a list of fields with the "Parent=" prefix

                checkIDs(line, feature_info, protein_id, "protein id parsing", globs);
                # A quick check to make sure we have read only one ID

                protein_id = protein_id[0].replace("protein_id=", "").replace("\"", "");
                # Unpack and parse the gene ID

                #globs['annotation'][gene_id]['transcripts'][parent_id]['protein-id'] = protein_id

                globs['tid-to-pid'][parent_id] = protein_id;

                # num_exons = len(globs['annotation'][gene_id]['transcripts'][parent_id]['exons']);
                # exon_id = "exon-" + str(num_exons+1);
                # # Because exon IDs are not always included for CDS, or they only represent the CDS as a whole (e.g. protein ID from Ensembl), we 
                # # count the number of exons in the transcript as the ID

                # globs['annotation'][gene_id]['transcripts'][parent_id]['exons'][exon_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand };
                # # Add the ID and related info to the annotation dict.

                num_prot_ids += 1;
                # Add to the number of exons read

        if num_prot_ids == 0:
            step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(num_prot_ids) + " protein IDs read");
        else:
            step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_prot_ids) + " protein IDs read");
        # Status update

        if len(missing_transcripts) == num_transcripts:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: " + str(len(missing_transcripts)) + " missing transcripts (present as parent of other feature, but not as its own feature)");      
        else:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(len(missing_transcripts)) + " missing transcripts (present as parent of other feature, but not as its own feature)");  
        
        if len(globs['tid-to-pid']) == 0:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: " + str(len(globs['tid-to-pid'])) + " protein IDs matched");
        else: 
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(len(globs['tid-to-pid'])) + " protein IDs matched");  
        #print("# INFO: " + str(len(missing_transcripts)) + " missing transcripts (present as parent of other feature, but not as its own feature)");
        # Report the number of missing transcripts

        # if num_cds_exons == 0:
        #     CORE.errorOut("GXF2", "No CDS exons found in input annotation file! Cannot calculate degeneracy without coding sequences.", globs);
        # # Check to make sure at least one CDS sequence is found, otherwise error out

    # Read coding exons to get protein IDs
    ####################

    elif globs['id-repl']:
        step = "Replacing strings in IDs according to -repl";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Status update

        ids_replaced = 0;
        for gene in globs['annotation']:
            for transcript in globs['annotation'][gene]['transcripts']:
                globs['tid-to-pid'][transcript] = transcript.replace(globs['id-repl'][0], globs['id-repl'][1]);
                # For every transcript in every gene, do the specified replacement and save it to the tid-to-pid list

                ids_replaced += 1;
                # Increment the counter

        if ids_replaced == 0:
            step_start_time = CORE.report_step(globs, step, step_start_time, "WARNING: " + str(ids_replaced) + " IDs edited");
        else:
            step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(ids_replaced) + " IDs edited");
        # Status update

    # Replace ID strings if -repl is specified
    ####################

    return globs;

    #############################################################################
    