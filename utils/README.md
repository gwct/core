# core/utils/

General purpose bash scripts.

| File | Description | 
| ------ | ----------- |
| backup | A script to quickly make a backup (`.bak`) of a file or folder by using `cp` and appending a timestamp to the backup. |
| cancel_all_jobs | This script looks up all currently running SLURM jobs for the current user and prompts to cancel them all at their own risk. Useful for when a large workflow needs to be cancelled and nothing else is running. |
| ensembl_get_fasta.sh | Takes a list of species and retrieves sequences from Ensembl. |
| fasize | This script calculates the size of one or more FASTA files in terms of the number of bases. It can also provide a more detailed analysis of the base composition and GC content. If a FASTA index file (.fai) is provided, it calculates the genome size based on the index (but cannot do further  counts of base composition). |
| git_add | Wraps `git add` and checks the file sizes of the provided files (or all changed files if none provided) to make sure none are over github's 100MB file size limit. If any are, they are listed and no files are staged. Prevents the situation when files over the limit are staged and committed, only to get the error on pushing that some are over the file size limit and having to un-commit and un-stage all files to deal with them |
| how_many_trees | A little script to show the number of possible rooted tree topologies for a given number of species. |
| interact | A bash command to start an interactive session on a SLURM partition. |
| jobres | This script wraps sacct to display job resource usage for a given user and/or job IDs. It parses the sacct output, which provides multiple records for each job ID, and calculates the maximum elapsed time and memory usage for each job.
| logcheck | Checks all files in a directory with a given extentsion (.log by default) for error and warning messages. |
| sres2 | This script displays the partition information for a SLURM cluster. |
| sresg | (SLOW, partitions hardcoded) This script displays the GPU information for particular partitions on a SLURM cluster. |
| sstart | Prints expected start times for all jobs of a given user stuck in a SLURM queue. |
