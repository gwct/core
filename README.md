# CORE: Code fOr Romps in Evolutionary data
### A mixture of scripts and libraries to help with sequence data manipulation, tree parsing, and other things.

## Author
#### Gregg Thomas

###### Please note that many of these scripts handle FASTA files. For my scripts, all FASTA files *must* have the extension .fa.
###### For any script, use the -h flag for specific usage details.

1. corelib/core.py
  * General helper functions such as reading sequences to a dictionary. You'll have to look to see what all is there.
2. corelib/treeparse.py
  * A couple functions that read newick formatted trees and return all relevant information in a more useful way to code with.
3. count_aln.py
  * This script gathers statistics about a single alignment file, or a directory full of alignment files.
4. count_pos.py
  * This script simply counts the number of amino acids or nucleotides in a file or directory.
5. fa_concat.py
  * Concatenates many FASTA formatted sequence files into a single FASTA file.
6. fa_edit.py
  * A general purpose FASTA handling script. Can relabel and trim headers and remove start and stop AAs.
7. how\_many\_trees
  * Just a little script to show the number of possible rooted tree topologies for a given number of species.
8. run_gblocks.py
  * A script to run GBlocks to mask a directory full of alignments in FASTA format. Note: This currently runs GBlocks at the most relaxed settings for phylogenetic tree inference. It will reject any masks that remove more than 20% of the columns from the original alignment.
9. run_muscle.py
  * This will make MUSCLE alignments out of a directory of FASTA files. Obviously, you'll need MUSCLE installed and in your PATH.
10. run_pasta.py	
  * This will make PASTA alignments out of a directory of FASTA files. Again, you'll need PASTA installed and in your PATH.
11. run_raxml.py
  * Runs some basic RAxML analyses on a directory full of FASTA files.
12. seq_convert.py
  * A sequence file format conversion tool. Currently converts between FASTA (.fa), Phylip (.ph), and Nexus (.nex) formats. It assumes files will have those extensions. Remember, these formats vary a lot in the details, so they might not work right away for everything. Let me know if you run into problems and I'll try to fix it.
