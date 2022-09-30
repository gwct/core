#############################################################################
# Functions for retrieving information from GTF and GFF files.
# Gregg Thomas
# Fall 2021-present
#############################################################################

import core
import gzip

#############################################################################

def checkIDs(l, info, id_list, step):
# Each time we read feature info and look for IDs, this checks to make sure only one
# ID is found. Probably unnecessary, but easy enough to check.

    if len(id_list) != 1:
        print("\n\n");
        print(id_list);
        print(info);
        print(l);
        print("\n\n");
        core.errorOut("GXF1", "Invalid number of IDs found during " + step);

#############################################################################

def getGenes(filename, file_type='gff', quiet=False):

    if not quiet:
        print("# Detecting compression of annotation file")
    compression = core.detectCompression(filename);
    if compression == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
        if not quiet:
            print("# Success: No compression detected");
    else:
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
        if not quiet:
            print("# Success: " + compression + " detected");
    # Detect the compression of the input annotation file

    if file_type == "gff":
        field_splitter = ";";
        gene_id_format = "ID=";
        transcript_id_format = "ID=";
        exon_id_format = "ID=";
        transcript_parent_format = "Parent=";
        exon_parent_format = "Parent=";

    elif file_type == "gtf":
        field_splitter = "; ";
        gene_id_format = "gene_id ";
        transcript_id_format = "transcript_id ";
        exon_id_format = "exon_id ";
        transcript_parent_format = "gene_id ";
        exon_parent_format = "transcript_id ";
    # These outline the differences between GFF and GTF

    annotation = {};
    # The main annotation storage dict. A nested structure for genes, transcripts, and coding exons.
    # <gene id> : { <header>, <start coord>, <end coord>, <strand>, 
    #                   { <transcript id> : <transcript header>, <transcript start coord>, <transcript end coord> <transcript strand>, 
    #                       { <exon id> : <exon header>, <exon start coord>, <exon end coord>, <exon strand> } } }

    ####################

    if not quiet:
        print("# Reading genes");

    with reader(filename) as open_file:
        for line in open_file:
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

                checkIDs(line, feature_info, feature_id, "gene id parsing");
                # A quick check to make sure we have read only one ID

                feature_id = feature_id[0].replace(gene_id_format, "").replace("\"", "");
                # Unpack and parse the ID

                annotation[feature_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand, 'transcripts' : {} };
                # Add the ID and related info to the annotation dict. This includes an empty dict for transcripts to be stored in a similar way
        ## End line loop
    ## Close file

    if not quiet:
        print("# Success: " + str(len(annotation)) + " genes read");
    # Status update

    return annotation, compression;

#############################################################################

def getTranscripts(filename, file_type="gff", quiet=False):

    annotation, compression = getGenes(filename, file_type, quiet);

    if compression == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
    else:
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
    # Get the appropriate reader fucntions based on compression

    if file_type == "gff":
        field_splitter = ";";
        gene_id_format = "ID=";
        transcript_id_format = "ID=";
        exon_id_format = "ID=";
        transcript_parent_format = "Parent=";
        exon_parent_format = "Parent=";

    elif file_type == "gtf":
        field_splitter = "; ";
        gene_id_format = "gene_id ";
        transcript_id_format = "transcript_id ";
        exon_id_format = "exon_id ";
        transcript_parent_format = "gene_id ";
        exon_parent_format = "transcript_id ";
    # These outline the differences between GFF and GTF

    if not quiet:
        print("# Reading transcripts");

    num_transcripts = 0;
    # A count of the number of transcripts read

    tid_to_gid = {};
    # A dictionary to link each transcript to it's gene ID, necessary for reading in exons in the next step
    # <transcript id> : <gene id>

    missing_genes = [];
    # Sometimes, transcripts will be present in the GXF file that have genes that don't have their own entry. This list
    # keeps track of those genes to report later

    with reader(filename) as open_file:
        for line in open_file:
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

                checkIDs(line, feature_info, feature_id, "transcript id parsing");
                # A quick check to make sure we have read only one ID

                feature_id = feature_id[0].replace(transcript_id_format, "").replace("\"", "");
                # Unpack and parse the ID

                parent_id = [ info_field for info_field in feature_info if info_field.startswith(transcript_parent_format) ];
                # Get the gene ID associated with the transcript as a list of fields with the "Parent=" prefix

                checkIDs(line, feature_info, parent_id, "transcript parent id parsing");
                # A quick check to make sure we have read only one ID

                parent_id = parent_id[0].replace(transcript_parent_format, "").replace("\"", "");
                # Unpack and parse the gene ID

                if parent_id not in annotation:
                    if parent_id not in missing_genes:
                        missing_genes.append(parent_id);
                    continue;
                # If the gene doesn't have it's own entry in the GXF file, add it to the list of missing genes and
                # skip this transcript
                else:
                    tid_to_gid[feature_id] = parent_id;
                # Otherwise, add the transcript id/gene id pair to the lookup dict

                annotation[parent_id]['transcripts'][feature_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand, 'exons' : {} };
                # Add the ID and related info to the annotation dict. This includes an empty dict for exons to be stored in a similar way

                num_transcripts += 1;
                # Add to the number of transcripts read
        ## End line loop
    ## Close file

    if not quiet:
        print("# Success: " + str(num_transcripts) + " transcripts read");
        # Status update

        print("# INFO: " + str(len(missing_genes)) + " missing genes (present as parent of other feature, but not as its own feature)");  

    return annotation, compression, tid_to_gid;

#############################################################################

def getExons(filename, file_type="gff", coding_only=False, quiet=False):

    annotation, compression, tid_to_gid = getTranscripts(filename, file_type, quiet);

    if compression == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
    else:
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
    # Get the appropriate reader fucntions based on compression

    if file_type == "gff":
        field_splitter = ";";
        gene_id_format = "ID=";
        transcript_id_format = "ID=";
        exon_id_format = "ID=";
        transcript_parent_format = "Parent=";
        exon_parent_format = "Parent=";

    elif file_type == "gtf":
        field_splitter = "; ";
        gene_id_format = "gene_id ";
        transcript_id_format = "transcript_id ";
        exon_id_format = "exon_id ";
        transcript_parent_format = "gene_id ";
        exon_parent_format = "transcript_id ";
    # These outline the differences between GFF and GTF

    if coding_only:
        if not quiet:
            print("# Reading coding exons");
        feature_str = "CDS";
    else:
        if not quiet:
            print("# Reading all exons");
        feature_str = "exon"

    num_cds_exons = 0;
    # A count of the number of transcripts read

    missing_transcripts = [];
    # Sometimes, exons will be present in the GXF file that have transcripts that don't have their own entry. This list
    # keeps track of those transcripts to report later

    pid_to_tid = {};
    # The ID of an exon might be linked to its transcript and this dictionary keeps track of that

    with reader(filename) as open_file:
        for line in open_file:
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

            if feature_type in [feature_str]:
                feature_id = [ info_field for info_field in feature_info if info_field.startswith(exon_id_format) ];
                # Get the feature ID as a list of fields with the "ID=" prefix

                #checkIDs(line, feature_info, feature_id, "exon id parsing", globs);
                # A quick check to make sure we have read only one ID

                if feature_id:
                    feature_id = feature_id[0].replace(exon_id_format, "").replace("\"", "");
                # Unpack and parse the ID

                parent_id = [ info_field for info_field in feature_info if info_field.startswith(exon_parent_format) ];
                # Get the transcript ID associated with the exon as a list of fields with the "Parent=" prefix

                checkIDs(line, feature_info, parent_id, "exon parent id parsing");
                # A quick check to make sure we have read only one ID

                parent_id = parent_id[0].replace(exon_parent_format, "").replace("\"", "");
                # Unpack and parse the gene ID

                if parent_id not in tid_to_gid:
                    if parent_id not in missing_transcripts:
                        missing_transcripts.append(parent_id);
                    continue;
                # If the transcript doesn't have it's own entry in the GXF file, add it to the list of missing transcripts and
                # skip this exon
                else:
                    gene_id = tid_to_gid[parent_id];
                # Otherwise, lookup the gene ID associated with the transcript to use below

                num_exons = len(annotation[gene_id]['transcripts'][parent_id]['exons']);
                exon_id = "exon-" + str(num_exons+1);
                # Because exon IDs are not always included for CDS, or they only represent the CDS as a whole (e.g. protein ID from Ensembl), we 
                # count the number of exons in the transcript as the ID

                pid_to_tid[feature_id] = parent_id;

                annotation[gene_id]['transcripts'][parent_id]['exons'][exon_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand, 'feature-id' : feature_id };
                # Add the ID and related info to the annotation dict.

                num_cds_exons += 1;
                # Add to the number of exons read
        ## End line loop
    ## Close file

    if not quiet:
        print("# Success: " + str(num_cds_exons) + " coding exons read");
        # Status update

        print("# INFO: " + str(len(missing_transcripts)) + " missing transcripts (present as parent of other feature, but not as its own feature)");  
        # Report the number of missing transcripts

    return annotation, compression, tid_to_gid, pid_to_tid;

#############################################################################