############################################################
# For substitution rates, 03.19
# Retrieves Ensembl orthologs for a given query species and
# a set of target species. Retrieves sequences.
############################################################

# homo_sapiens
# aotus_nancymaae,bos_taurus,chlorocebus_sabaeus,gorilla_gorilla,macaca_mulatta,microcebus_murinus,monodelphis_domestica,mus_musculus,nomascus_leucogenys,pan_troglodytes,pongo_abelii
# Human

# aotus_nancymaae
# homo_sapiens,bos_taurus,chlorocebus_sabaeus,gorilla_gorilla,macaca_mulatta,microcebus_murinus,monodelphis_domestica,mus_musculus,nomascus_leucogenys,pan_troglodytes,pongo_abelii


import core, sys, os, argparse, treeparse as tp, requests
import xml.etree.ElementTree as ET

print "Program call: " + " ".join(sys.argv) + "\n"

parser = argparse.ArgumentParser(description="");
parser.add_argument("-q", dest="query", help="The path to a FASTA file with sequences from the target species.", default=False);
parser.add_argument("-n", dest="query_name", help="The name of the target species.", default=False);
parser.add_argument("-t", dest="targets", help="A comma separated list of species to query.", default=False);
parser.add_argument("-m", dest="mode", help="Which orthologs to retrieve: 'all' (default), 'oto'.", default='all');
parser.add_argument("-s", dest="seqtype", help="Which type of sequence to retrieve: 'cds' (default), 'prot'.", default='cds');
parser.add_argument("-o", dest="output", help="The output file name.", default=False);
args = parser.parse_args();
# Input options.

if not any(opt for opt in [args.query, args.query_name, args.targets, args.output]):
    sys.exit(" * ERROR: All of -q, -n, -t, and -o must be specified.");
if not os.path.isfile(args.query):
    sys.exit(" * ERROR: Cannot find specified query fasta file (-q): " + args.query);
if args.mode not in ['all', 'oto']:
    sys.exit(" * ERROR: -m must be either 'all' or 'oto'.");
if args.seqtype not in ['cds', 'prot']:
    sys.exit(" * ERROR: -s must be either 'cds' or 'prot'.");
targets = args.targets.split(",");
##############

start = core.getLogTime();
logfilename = args.output.replace(".txt","") + "-" + start + ".log";

server = "http://rest.ensembl.org"
# The Ensembl REST server.

ordered_speclist = args.query_name.split() + targets

with open(args.output, "w") as outfile, open(logfilename, "w") as logfile:
    logfile.write(start + "\n");
    logfile.write("Program call: " + " ".join(sys.argv) + "\n\n");

    outfile.write("\t".join(ordered_speclist) + "\n");
    query_seqs = core.fastaGetDict(args.query);
    num_orths = 0;

    numbars, donepercent, i, numseq = 0, [], 0, len(query_seqs);
    for title in query_seqs:
        numbars, donepercent = core.loadingBar(i, numseq, donepercent, numbars);
        i += 1;
        # Loading bar

        orths = { n : [] for n in ordered_speclist };
        # Initialize the current orth dictionary.

        gid = title[title.index("gene:"):]
        gid = gid[gid.index(":")+1:gid.index(".")];
        # Get the gene ID to query.

        #tid = title[title.index("_")+1:title.index(" ")];
        #query_seq = query_seqs[title];

        orths[args.query_name].append(gid);
        # Add the query gene ID into the orth dict.

        #gid = "ENSG00000188086";
        #print gid;
        logfile.write(gid + "\n");

        ext = "/homology/id/" + gid + "?type=orthologues;sequence=cds;aligned=0;format=condensed;"
        # Make the Ensembl query

        for t in targets:
            ext = ext + "target_species=" + t + ";";
        ext = ext[:-1];
        # Add the target species to the query

        attempt = 2;
        r = requests.get(server+ext, headers={ "Content-Type" : "text/xml"})
        while not r.ok:
            logfile.write(" -> Attempt " + str(attempt) + "\n");
            r = requests.get(server+ext, headers={ "Content-Type" : "text/xml"});
            if attempt > 100:
                #print r.raise_for_status()
                print " -> Coult not connect... skipping...";
                continue;
                #sys.exit()
        # Query as XML and retrieve from Ensembl

        root = ET.fromstring(r.content);
        #print root.tag;
        #print r.text;
        # Parse XML

        for homology in root.iter('homologies'):
            #print homology.tag, homology.attrib;
            orths[homology.attrib['species']].append(homology.attrib['id']);
        # Get all the gene IDs for every homologous gene.

            # for sub in homology:
            #     if sub.tag == 'source':
            #         orths[args.query_name][0]['protein'] = sub.attrib['protein_id'];
            #         if sub.attrib['seq'] != query_seqs[title]:
            #             print " * Warning: query seqs not identical for: " + gid;
                        
            #             sys.exit();
            #     # Get the protein ID from the query and ensure that the query seq from the file and from the database are the same.

            #     if sub.tag == 'target':
            #         print " -->", sub.attrib['species'];
            #         orths[sub.attrib['species']].append( { 'gene' : sub.attrib['id'], 'protein' : sub.attrib['protein_id'], 'seq' : sub.attrib['seq'] } );
            #     # Get the target
            # This chunk was for when I thought I could get sequences at the same time, but it only gets 1 transcript, often not the one I need (longest).

        #print orths;
        #print "----------";

        outline = "";
        if args.mode == 'oto' and all(len(orths[spec]) == 1 for spec in orths):
            num_orths += 1;
            for spec in ordered_speclist:
                outline += orths[spec][0] + "\t";
            outfile.write(outline[:-1] + "\n");

        else:
            num_orths += 1;
            for spec in ordered_speclist:
                if orths[spec] != []:
                    outline += ",".join(orths[spec]) + "\t";
                else:
                    outline += "\t";
            outfile.write(outline[:-1] + "\n");

        #print "------------\n\n";
        #sys.exit();
    logfile.write("----------\n");
    logfile.write(str(num_orths) + " otos written\n");

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print num_orths, "orths written.";
