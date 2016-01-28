
# baits-analysis
scripts for getting full paralog sequences from fastq files from bait sequencing

The final plan is for the scripts to be in two parts.  One of them will be accessed through tfastq_assembly_master.py and the other through tcontig_classif_master.py.  The second of these scripts is finished, but the first is still in progress.  The scripts that it will eventually use are all finished and can be called individually, though.

trans_bcparse_2reads.py
This parses the original fastq files according to their inline barcodes into separate files for each individual.
It is not part of either of the master scripts.

##############################################
tfastq_assembly_master.py--It is not yet written, although the scripts it calls are.
The master script that takes the fastq files that are separated according to individual and makes contigs that are split according to gene family.
It calls the following scripts, in addition to blast and the assembly program of choice (spades, mazurca, others?)

trans_fastq_to_2_blast.py
This makes separate fasta files from the fastq files and writes the blast script.

tbaits_blastn_parse.py
This parses the blast results and to determine which sequences belong to which gene family.

tblast_to_fastq.py
This makes separate fastq files for each gene family, each of which contains sequences from multiple individuals.
It writes a script for spades to assemble the data.

tassembly_to_blast.py
This parses the spades output and makes a new blast database that includes the old sequences and the new spades contigs.
It also writes the script to blast the original fasta files against this database.

tbaits_blastn_parse.py
(again)--something is missing for spades here***

tassembly_to_loci.py
This parses the spades output and makes separate fasta files from the contigs to blast them against the original transcriptome sequences
to determine which are the exons and which are the introns.
It also writes the blast scripts.

############################################
tcontig_classif_master.py
The master script that takes the output from blasting the spades contigs against the original transcriptomes, separates the exons,
and assembles them into full sequences for each paralog from each gene family.
It calls the following scripts, in addition to mafft and raxml.

fasta_to_phylip.py
This script converts the fasta file produced by mafft to a phylip file that is readable by raxml.  It is used multiple times.

tbaits_intron_removal.py
This removes the introns from the contigs based on results of blasting them against transcriptome sequences.
It also writes a script to use mafft to align the contigs to the backbone alignments and add the short sequences to raxml trees.

tseq_placer_dup.py
This groups the contigs according to their position in the raxml trees and labels each group according to paralog.

tcontigs_to_fixed_paralogs.py
This tried to combine the sequences from each individual in each group into one full sequence for that paralog.
It writes a script to align the unclassifiable contigs to the backbone alignment and add them to the backbone trees for the next 
round of analysis.

tparalog_combiner.py
This combines the new sequences from the various clades together with the backbone sequences to make new backbone trees.
It also makes various information files.

tseq_placer_dup.py, tcontigs_to_fixed_paralogs.py, and tparalog_combiner.py are repeated for several rounds.

At the end of this, the following scripts are run:

(tbaits_cleanup.py)--not added yet

tundivcontigs_combiner.py
This takes all of the contigs for each locus that could not be successfully merged into a complete sequence and combines them into one file per locus.

tcontig_selection.py
This looks through the files for each locus and selects the best contigs to add them to the tree.

tparcomb_combiner.py
This script combines all of the sequence and information files from the previous rounds.  It also looks for individuals that have more than one sequence per
paralog, when those sequences were found in different rounds.  It flags these sequences.

tparcomb_final.py
This script attempts to combine the sequences flagged by tparcomb_combiner.py into a single sequence.  It writes a script to make alignments and trees
with the final set of sequences.  It also makes a list of ambiguous sequences that could not be classified according to paralog.

tambig_seq_renaming.py
After looking at the trees by hand and classifying the ambiguous sequences that could be classified, this script renames those sequences in the
alignment files, the tree files, and the dictionary that says which sequences belong to which paralogs.
