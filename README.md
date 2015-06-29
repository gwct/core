# core
### A mixture of scripts and libraries to help with sequence data manipulation, tree parsing, and other things.

## Author
#### Gregg Thomas

###### Please note that many of these scripts handle FASTA files. For my scripts, all FASTA files *must* have the extension .fa.
###### For any script, use the -h flag for specific usage details.

1. concatenate_fa.py
  * Concatenates many FASTA formatted sequence files into a single FASTA file.
2. corelib/core.py(c)
  * General helper functions such as reading sequences to a dictionary. You'll have to look to see what all is there.
3. corelib/treeparse.py(c)
  * A couple functions that read newick formatted trees and return all relevant information in a more useful way to code with.
4. count_aln.py
  * This script gathers statistics about a single alignment file, or a directory full of alignment files.
5. count_pos.py
  * This script simply counts the number of amino acids or nucleotides in a file or directory.
6. fa2nex.py
  * A FASTA to Nexus converter.
7. fa2phy.py
  * A FASTA to Phylip converter.
8. fa_edit.py
  * A general purpose FASTA handling script. Can relabel and trim headers and remove start and stop AAs.
9. how\_many\_trees
  * Just a little script to show the number of possible rooted tree topologies for a given number of species.
10. nj_tree.r
  * This will create a neighbor joining tree from a given distance matrix in a certain format. I made this for a very specific purpose, so it's probably not too useful at this point...
  * Usage:
  * `rscript nj_tree.r [input file with distance matrix] [outgroup for rooting] [output file name]`
11. phy2fa.py
  * A Phylip to FASTA converter.
12. run_muscle.py
  * This will make MUSCLE alignments out of a directory of FASTA files. Obviously, you'll need MUSCLE installed and in your PATH.
13. run_pasta.py	
  * This will make PASTA alignments out of a directory of FASTA files. Again, you'll need PASTA installed and in your PATH.
