# core/python/

Python scripts for many tasks including sequence handling, tree making, and sequence alignment.
Some of these programs are mainly used as wrappers to easily run other genomics or phylogenetics programs on a bunch of files. 

All of these scripts are written in Python 3.4+ (https://www.python.org/downloads/).

### For any script, use the -h flag for specific usage details.

| Directory/File | Description | 
| ------ | ----------- |
| generators/ | Scripts that generate commands to run on many files as well as SLURM submission scripts. Includes generators for a coding sequence alignment pipeline. |
| lib/| Functions for many basic tasks as well as specific sequence or tree tasks. |
| python-template/ | A folder that includes a template for a standalone Python-based project. |
| ensembl_get_orths.py | This script uses Ensembl's Rest API to get ortholous genes from a set of species. |
| fastand_furious.py | A general purpose FASTA handling script. Can count positions in sequences and alignments, concatentate alignments, combine and split .fa files, trim and relabel headers, remove sequences and start positions, and replace specific states in sequences. |
| fastq_stats.py | Summarizes reads in FASTQ files. |
| gxf_feature_counter.py | Counts features in GTF or GFF files. Needs further development. |
| gxf_parse.py | Converts certain features in GTF or GFF files into a more bed-like tab-delimited format. |
| seq_convert.py | Can convert sequences between FASTA (.fa), Phylip (.ph), and Nexus (.nex) formats. Note that these formats often vary in small ways between users, so this might not work right away for you. Consider this in Beta. |
| tree.py | Some general purpose Newick tree handling modules. Can join or separate directories or files of trees, label internal nodes of trees, check if trees are rooted, root trees, calculate concordance factors, and count and relabel tips. Dependencies: Newick Utilities (http://cegg.unige.ch/newick_utils), called as `nw_reroot` is required to root trees with `--root`. |