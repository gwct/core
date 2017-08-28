# CORE: COde for Romps in Evolutionary data
### A mixture of scripts and libraries to help with sequence data manipulation, tree parsing, and other things.

## Author
#### Gregg Thomas

## About
###### These scripts can be used for many tasks including sequence handling, tree making, and sequence alignment.
###### Some of these programs are mainly used as wrappers to easily run other genomics or phylogenetics programs on a bunch of files. Pay attention to the dependencies for each script to make sure you have the proper programs installed.
###### Please note that many of these scripts expect input as FASTA files. For my scripts, these *must have the extension .fa*. If you don't have FASTA formatted files, you can use seq_convert to get them to FASTA format and fa_edit to make any changes you need to them afterwards.

###### Almost all of these scripts are written in Python 2.7 (https://www.python.org/downloads/).
###### For any script, use the -h flag for specific usage details.

## CORE scripts

1. corelib/core.py
  * General helper functions such as reading sequences to a dictionary. You'll have to look to see what all is there.
2. corelib/fastalib.py
  * All the modules for the fasta_nd_furious.py script for FASTA handling are here.
3. corelib/nj_tree.r
  * Simple R script to get a Neighbor Joining tree. Used by wrappers --sdm and probably not helpful standalone.
  * Dependencies:</br>
&emsp;i. R (https://www.r-project.org/)
4. corelib/treeparse.py
  * Lots of functions for Newick trees. Used by a lot of my scripts.
5. corelib/wrapperlib.py
  * All the modules to run the other programs in wrappers.py
6. dev/
  * Some scripts that I use, but only occasionally.
7. legacy-scripts/
  * The old versions of FASTA handling scripts (before fasta_nd_furious.) and wrappers for other programs (before wrappers.py).
8. fasta_nd_furious.py
  * A general purpose FASTA handling script. Can count positions in sequences and alignments, concatentate alignments, combine and split .fa files, trim and relabel headers, remove sequences and start positions, and replace specific states in sequences.
9. flight_of_the_concordance.py
  * Given a set of single-copy gene trees and a species tree, this script calculates concordance factors for each node in the species tree.
10. how\_many\_trees
  * Just a little script to show the number of possible rooted tree topologies for a given number of species.
11. isofilter.py
  * A script for filtering all but the longest isoform from a set of proteins. Works with Ensembl and NCBI data.
12. paml_lrt.py
  * After using wrappers.py to run the branch-site test with codeml, this script does the likelihood ratio tests for all genes and reports those that pass.
13. root_trees.py
  * A simple script to re-root a set of gene trees given an outgroup species.
  * Dependencies:</br>
&emsp;i. Newick Utilities (http://cegg.unige.ch/newick_utils), called as `nw_reroot`</br>
14. seq_convert.py
  * Can convert sequences betwee FASTA (.fa), Phylip (.ph), and Nexus (.nex) formats.
15. wrappers.py
  * This script is big. It is meant as a wrapper for all of the other evolutionary analysis programs I use, including MUSCLE, PASTA, GBlocks, RAxML, codeml, SDM, and r8s. Many of these programs are meant to work on a single file, or require an input control file, or the data to be rearranged between steps. This script hopefully does that for all of these programs (for the analyses I usually do). Obviously, to run one of the modules you will need that program installed!

##