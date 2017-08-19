
# baits-analysis
scripts for getting full paralog sequences from fastq files from bait sequencing

The scripts are in three major parts: the first accessed through tfastq_assembly_master.py, the second accessed through tcontig_classif_master.py, and the third accessed through tgenefam_to_spptree_master.py.  The scripts that these master scripts call can also be run individually.

All scripts have documentation and trying to run one of the scripts without additional arguments/with the wrong number of arguments will give you the list of arguments for that script.

All master scripts use gnu-parallel and they all write files for the slurm workload management system.  If they are run in Array mode (when that option is available), they can only be run on slurm systems.  
If they are run in parallel mode, they should be able to run on other workload management systems, although this has not been tested.

The scripts are written for python 2.7.  Required python modules include:
Bio
collections (defaultdict, Counter, OrderedDict)
copy
dendropy (If you are getting errors where dendropy does not recognize some of the commands, you probably need to update to the most recent version of dendropy)
gzip (If you are getting errors reading gzipped files directly, gunzip the files yourself first and it will work with the gunzipped files.)
itertools (izip)
numpy
random
re
subprocess
sys

Other required programs include:
blast
gnu-parallel
mafft
notung
raxml
spades

Optional programs include:
mrbayes
minimo
ssake

trans_bcparse_2reads.py
This parses the original fastq files according to their inline barcodes into separate files for each individual.
It is not part of any of the master scripts.

##############################################
tfastq_assembly_master.py
The master script that takes the fastq files that are separated according to individual and makes contigs that are split according to gene family.
It calls the following scripts, in addition to blast and the spades assembler.
The second assembly step can also be performed with the ssake or minimo assemblers.
The output of this master script is a pair of files for each gene family for each group of organisms: a fasta file with the assembled contigs (labeled according to individual) and a blast output file showing how those contigs blasted to the sequence database (and thus showing which parts of those sequences are introns and which are exons, if the database does not contain intron data). 

trans_fastq_to_2_blast.py
This makes separate fasta files from the fastq files and writes the blast script.

##Round 1 (optional, but should be run if there are few input sequences per individual or those sequences are quite divergent from the sequences in the blast database):

**blast search**

tbaits_blastn_parse.py
This parses the blast results to determine which sequences belong to which gene family.

tblast_to_fastq.py (together mode)
This makes separate fastq files for each gene family, each of which contains multiple individuals, and writes a script for spades to assemble the reads.

**spades assembly**

tassembly_to_blast.py (spades mode)
This parses the spades output and makes a new blast database that includes the old sequences and the new spades contigs.
It also writes the script to blast the original fasta files against this database.

##Round 2 (or the only round, if it is run for a single round):

**blast search**

tbaits_blastn_parse.py
This parses the blast results to determine which sequences belong to which gene family.

tblast_to_fastq.py (separate mode)
This makes separate fastq files for each gene family, with a separate fastq file for each individual/gene family combination, and writes a script for spades to assemble the reads.

**spades assembly**

tassembly_to_loci.py (spades mode)
This parses the spades output and makes separate fasta files from the contigs to blast them against the original transcriptome sequences to determine which are the exons and which are the introns.
It also writes the blast scripts.

**blast search**

Separate scripts are also written for cleanup and to run Round 2 using the minimo or ssake assemblers, instead of spades.

##Round 2 with minimo:

tblast_to_fasta_qual.py (separate mode)
This makes separate fasta files and base quality files for each gene family, with a separate pair of files for each individual/gene family combination, and writes a script for minimo to assemble the reads.

**minimo assembly**

tassembly_to_loci.py (minimo mode)
This parses the minimo output and makes separate fasta files from the contigs to blast them against the original transcriptome sequences to determine which are the exons and which are the introns.
It also writes the blast scripts.

**blast search**

##Round 2 with ssake:

tblast_to_ssake.py (separate mode)
This makes separate fastq files for each gene family, with a separate fastq file for each individual/gene family combination, and writes a script for ssake to convert the reads to the proper format and assemble them.

**ssake assembly**

tassembly_to_loci.py (ssake mode)
This parses the ssake output and makes separate fasta files from the contigs to blast them against the original transcriptome sequences to determine which are the exons and which are the introns.
It also writes the blast scripts.

**blast search**

############################################
tcontig_classif_master.py
The master script that takes the output from blasting the spades contigs against the original transcriptomes, separates the exons,
and assembles them into full sequences for each paralog from each gene family.
It calls the following scripts, in addition to mafft and raxml.

##Used multiple times:

fasta_to_phylip.py
This script converts the fasta file produced by mafft to a phylip file that is readable by raxml.

##Before the rounds of sequence classification begin:

tbaits_intron_removal.py
This removes the introns from the contigs based on results of blasting them against transcriptome sequences.
It also writes a script to use mafft to align the contigs to the backbone alignments and use raxml to place the contigs in the backbone trees

##The following are run for multiple rounds of sequence classification (generally 6, although the number can be varied):

If tcontig_classif_master.py is run in Array mode, these scripts are written by the tccm_group_array_subscript.py script.
If tcontig_classif_master.py is run in Parallel mode, it writes these scripts directly.

**mafft alignment**

**raxml sequence placement**

tseq_placer_dup.py
This groups the contigs according to their position in the raxml trees and labels each group according to paralog.

tcontigs_to_fixed_paralogs.py
This tried to combine the sequences from each individual in each group into one full sequence for that paralog.
It writes a script to align the unclassifiable contigs to the backbone alignment and add them to the backbone trees for the next 
round of analysis.

tparalog_combiner.py
This combines the new sequences from the various clades together with the backbone sequences to make new backbone trees.
It also makes various information files.

**raxml to create the new backbone gene trees** (skipped in the last round)

##Making the final alignments after the last round of sequence classification:

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

############################################
tgenefam_to_spptree_master.py
This master script takes the sequences that were created by tcontig_classif_master.py, makes gene family trees with them, and breaks those trees up into individual sets of sequences for each locus for phylogenetic analysis.
This script requires a species tree.  We have generally reconstructed the original species tree from chloroplast and nrDNA genes that were captured incidentally during the targeted sequence capture process.

tcombpars_to_trees.py
This script takes the sets of sequences from each gene family and writes a script to align them using mafft, build gene family trees using raxml, and use notung to find duplications in those trees.

**mafft alignment**

**raxml to create gene family trees**

##Round 1

**notung to find duplications, using the species tree provided**

tnotung_homolog_parsing.py
This script parses the notung results to determine which duplications are likely real, makes separate sequence files for the sequences from each paralog, and writes a script to align those sequences.

**mafft alignment**

tparalog_selector.py
This script analyzes the output from tnotung_homolog_parsing.py to make lists of loci and individuals that should be included in a final set of alignments, depending on the desired amount of missing data.

tal_combiner.py
This script takes a list of loci and a list of individuals and makes sequence alignments in various formats and the script to make a species tree (or the trees to be used in creation of that species tree) with those alignments.
Four modes are currently available:
combined: concatenated alignment for all loci and a raxml tree
combinednobs: concatenated alignment for all loci and a raxml tree without the bootstrap search (support values are not needed for the species tree for the 2nd round)
separate: separate alignments and raxml trees for each locus
separatemb: separate alignments and MrBayes trees for each locus (using fasta_to_nexus_mb.py)

At the end of Round 1, that species tree can be accepted (if the original species tree was good; Good mode for tgenefam_to_spptree_master.py), the concatenated tree from raxml can be used as the species tree for the next round (New mode for tgenefam_to_spptree_master.py), or the script can stop after tal_combiner.py runs to allow the user to create a species tree using their program of choice as the species tree for the next round (Pause mode for tgenefam_to_spptree_master.py).


###########################################
Scripts for making the original species tree from a set of single copy loci that are also present in the sequence reads (e.g., chloroplast, mitochondrial, or nrDNA sequences):

torig_spp_tree_master.py
This script is similar to tfastq_assembly_master.py in that it runs a separate blast search for the sequences used to build the species trees, if they were not in the original blast database.
It also uses blast and spades, along with the following scripts that are also used by tfastq_assembly_master.py: trans_fastq_to_2blast.py, tbaits_blastn_parse.py, tblast_to_fastq.py, and tassembly_to_loci.py
Only the first and the last scripts are different:

torig_spp_tree_blasting.py
This writes the script to blast the sequences against the blast database.

tbaits_to_spptreeseqs.py
After the spades assembly has been blasted back to the species tree sequence blast database, this script chooses the longest contig of the ones that blasted to each gene as the one to use for that individual for that gene.


After the master script has been run to find the sequences for each group, they will need to be combined into a single file for each gene.  This can be done using
tcontig_moving.py (Seqs mode).

The files of sequences will then need to be aligned.  After they are aligned, tal_combiner.py can be used to build separate or combined trees for the various loci.


###########################################
Other potentially useful scripts:

##organizing/renaming alignments and trees

alignment_combiner.py
This script combines two alignments while keeping only one copy of duplicate sequences.  It can input or output any format that is recognized by Biopython.

fasta_renaming.py
This script specifically renames fasta files.  If there is a long list of names that need to be changed, it runs much faster than treerenamer.py.

seqsorting.py
This script sorts the sequences in an alignment alphabetically.  It can input or output any format that is recognized by Biopython.

treerenamer.py
Given a list of current names and their replacement values, this script makes the replacements in any text file.

##figure making

tdup_mapping.py
This script uses the output from tnotung_homolog_parsing.py to map the gene duplications onto a species tree as branch labels.

tfig_making_prep.py
This script uses the output from tnotung_homolog_parsing.py and tparalog_selector.py and formats it so that R can use it for making coverage figures.

tfig_making.r
This script makes various coverage figures using the output of tfig_making_prep.py

heatmap_making.r
This script makes heatmaps of coverage per locus across individuals.
