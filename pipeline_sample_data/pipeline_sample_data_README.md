Files:

original sequences:
in the Anacampserotaceae, Basellaceae, and Cactaceae folders: two fastq files per individual: t1_[individual_name_number]_R1.fastq and t1_[individual_name_number]_R2.fastq 

for tfastq_assembly_master.py:

blast database: sample_blastdb.*
individual blast databases for the 3 gene families, for finding exons: exon_dbs (folder)
list of other names for various loci: sample_loci_shortened.txt
list of individuals for the two groups: Anac_list.txt, Bas_list.txt, Cact_list.txt (in the Anacampserotaceae, Basellaceae, and Cactaceae folders, respectively)


for tcontig_classif_master.py:

list of loci to be used in the analysis: sample_locus_list.txt
list of other names for various loci: sample_loci_shortened.txt
list of paralogs that belong to the various loci: sample_pardict.txt
list of outgroups to root the gene family trees: sample_outgroup_list.txt
list of the full group names and their abbreviations: sample_group_list.txt
no polyploids, so "none" instead of PolyList
list of individuals together with their groups: sample_ind_list_groups.txt
folder with the backbone trees: backbone_trees/
prefix for the backbone trees and alignments: new_
suffix for the backbone alignments: _al


for tgenefam_to_spptree_master.py:

list of loci to be used in the analysis: sample_locus_list.txt
list of individuals together with their groups, this time including the outgroups: sample_ind_list_groups_og.txt
species tree: RAxML_bestTree.sample_combined_22inds_4seqs

#Note that the master scripts do not start anything on their own.  They write scripts that you have to then start by hand, either by typing the script path, if you are not
#running on slurm or by typing sbatch /script/path if you are running them on slurm.

#Sample scripts for running with slurm:

***The most important command if you are running this on slurm!!***

egrep "Killed|Error|Check|CANCELLED|error|ERROR|exiting" slurm-*

There should never be a failed Dependency or a python error.  If you get either of these, even if it seems as if everything ran correctly, it did not, and you need to find the problem.
Once you have checked for errors (assuming there are no more running scripts), you can delete the slurm output files before starting the next part of the pipeline.
Part Three expects to have the Notung-2.8.1.6-beta.jar script in the scripts folder as well.  The remaining programs are presumed to be on the server and accessible normally.  If you get an error message associated with any of the programs, it is likely not present or called something else on your cluster/computer.

**Part One**
/users/ajm3/data/ajm3/scripts/tfastq_assembly_master.py 1June2017 /users/ajm3/data/ajm3/scripts/ /users/ajm3/data/ajm3/pipeline_sample_data/Anacampserotaceae/ t1_ /gpfs/scratch/ajm3/slurm/Anacampserotaceae/ Anac /users/ajm3/data/ajm3/pipeline_sample_data/exon_dbs/ none /users/ajm3/data/ajm3/pipeline_sample_data/Anacampserotaceae/Anac_list.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_loci_shortened.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_blastdb 4 Array 2
/users/ajm3/data/ajm3/scripts/tfastq_assembly_master.py 1June2017 /users/ajm3/data/ajm3/scripts/ /users/ajm3/data/ajm3/pipeline_sample_data/Basellaceae/ t1_ /gpfs/scratch/ajm3/slurm/Basellaceae/ Bas /users/ajm3/data/ajm3/pipeline_sample_data/exon_dbs/ none /users/ajm3/data/ajm3/pipeline_sample_data/Basellaceae/Bas_list.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_loci_shortened.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_blastdb 4 Array 2
/users/ajm3/data/ajm3/scripts/tfastq_assembly_master.py 1June2017 /users/ajm3/data/ajm3/scripts/ /users/ajm3/data/ajm3/pipeline_sample_data/Cactaceae/ t1_ /gpfs/scratch/ajm3/slurm/Cactaceae/ Cact /users/ajm3/data/ajm3/pipeline_sample_data/exon_dbs/ none /users/ajm3/data/ajm3/pipeline_sample_data/Cactaceae/Cact_list.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_loci_shortened.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_blastdb 4 Array 2

The following are what you will need to change to run it on your cluster:
/users/ajm3/data/ajm3/scripts/ is the folder in which all of the scripts are found.
1June2017 is the date.
/users/ajm3/data/ajm3/pipeline_sample_data/ is the folder in which all of the sample data are found
/users/ajm3/data/ajm3/pipeline_sample_data/Anacampserotaceae/ /users/ajm3/data/ajm3/pipeline_sample_data/Basellaceae/ and /users/ajm3/data/ajm3/pipeline_sample_data/Cactaceae/ are the folders in which the fastq files are found.
/gpfs/scratch/ajm3/slurm/Anacampserotaceae/ /gpfs/scratch/ajm3/slurm/Basellaceae/ and /gpfs/scratch/ajm3/slurm/Cactaceae/ are the scratch folders in which the analyses are actually run.  **Note that you must make this folder; it does not try to make it itself.**

The script to run (for me; in the case of Anacampserotaceae) is /gpfs/scratch/ajm3/slurm/Anacampserotaceae/fastq_assembly1.sh

When you check for errors, you will find things with the word "error" in them with spades.  But as long as you don't have any python errors, you're fine.
For this reason, you can use the alternative for this round **only**.:
egrep "Killed|Check|CANCELLED|ERROR|exiting" slurm-* 

**Part Two**

/users/ajm3/data/ajm3/scripts/tcontig_classif_master.py 6June2017 /users/ajm3/data/ajm3/scripts/ /gpfs/scratch/ajm3/slurm/ s2_final_ spades_ s /gpfs/scratch/ajm3/slurm/output/ slurm /users/ajm3/data/ajm3/pipeline_sample_data/sample_locus_list.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_loci_shortened.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_pardict.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_outgroup_list.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_group_list.txt none /users/ajm3/data/ajm3/pipeline_sample_data/sample_ind_list_groups.txt /users/ajm3/data/ajm3/pipeline_sample_data/backbone_trees/ new_ _al 8 6 Array

The following are what you will need to change to run it on your cluster:
/users/ajm3/data/ajm3/scripts/ is the folder in which all of the scripts are found.
6June2017 is the date.
/users/ajm3/data/ajm3/pipeline_sample_data/ is the folder in which all of the sample data are found (including the fastq files).
/gpfs/scratch/ajm3/slurm/ is the folder in which the folders with the sequences (e.g., Anacampserotaceae) are found.
/gpfs/scratch/ajm3/slurm/output/ is the folder where the output should go; you must also make this folder first.
slurm is the prefix for the output.  The fact that it is called slurm is purely coincidental.  The thing that makes this script run for slurm is the Array at the end.

The script to run (for me) is /gpfs/scratch/ajm3/slurm/output/slurm_contig_classification_script1.sh

**Part Three**

/users/ajm3/data/ajm3/scripts/tgenefam_to_spptree_master.py 6June2017 /users/ajm3/data/ajm3/scripts/ /gpfs/scratch/ajm3/slurm/output/ slurm _allover150_al.fa _o150 /users/ajm3/data/ajm3/pipeline_sample_data/sample_locus_list.txt /users/ajm3/data/ajm3/pipeline_sample_data/sample_ind_list_groups_og.txt 4 /users/ajm3/data/ajm3/pipeline_sample_data/RAxML_bestTree.sample_combined_22inds_4seqs Oryza_sativa_cds Array New

The following are what you will need to change to run it on your cluster:
/users/ajm3/data/ajm3/scripts/ is the folder in which all of the scripts are found.
1June2017 is the date.
/users/ajm3/data/ajm3/pipeline_sample_data/ is the folder in which all of the sample data are found (including the fastq files).
/gpfs/scratch/ajm3/slurm/output/ is the folder where the output should go.  It will make this folder.
slurm is the prefix for the output, still purely coincidentally.
Part Three expects to have the Notung-2.8.1.6-beta.jar script in the scripts folder as well.  The remaining programs are presumed to be on the server and accessible normally.  If you get an error message associated with any of the programs, it is likely not present or called something else on your cluster/computer.
If you have a different version of Notung, you will need to change the version number in both tcombpars_to_trees.py for round 1 and in tgenefam_to_spptree_master.py for round 2).

The script to run (for me) is /gpfs/scratch/ajm3/slurm/output/slurm_Species_Tree_Making1.sh

When tgenefam_to_spptree_master.py is run in "New" mode (as here), it will run for two rounds, using the given species tree (from ITS, matK, ndhF, and rbcL sequences) in the first round and a tree made from concatenated sequences of the new loci in the second round.

To run Astral (from within the /gpfs/scratch/ajm3/slurm/output/slurm_prunedgenetrees2/slurmpgts2_raxmlbs/ folder (changing the version of astral to the current one and the path to the correct path for astral):
java -jar ~/source/ASTRAL-master/Astral/astral.4.10.10.jar -i slurmpgts2_bestTrees.tre -b bootstrap_filelist.txt -r 100 -o slurmpgts2_astral_100bs.tre
The tree with bootstrap values and (internal) branch lengths will be the last one in the resulting slurmpgts2_astral_100bs.tre file.

#Sample scripts for running with bash commands.  (These scripts can also be run in slurm, because they have the proper slurm headers, but can be run on a desktop or (likely with some modification) on a non-slurm cluster.)

The way these scripts are written assumes that the python scripts are all in your path.  If this is not true (or you do not know what this means), and you need to provide the path to the scripts, replace the first "none" with the path to the folder the scripts are in and write the full path to the three master scripts at the start (i.e., /home/abby/scripts/tfastq_assembly_master.py instead of just tfastq_assembly_master.py).

When not run on a slurm cluster, these scripts will give you error messages such as:
/home/abby/transcriptomes/shell/Anacampserotaceae/7June2017_fastq_assembly_master.sh: line 7: module: command not found
These can be ignored, as these programs will be available all of the time, without needing to load them, on your personal computer.

In addition, all messages about "mkdir: cannot create directory" or "rm: file not found" can also be ignored.

Any RAxML or python error messages indicate that something is wrong with the pipeline and must be dealt with.

**Part One**
tfastq_assembly_master.py 7June2017 none ~/transcriptomes/pipeline_sample_data/Anacampserotaceae/ t1_ ~/transcriptomes/shell/Anacampserotaceae/ Anac ~/transcriptomes/pipeline_sample_data/exon_dbs/ none ~/transcriptomes/pipeline_sample_data/Anacampserotaceae/Anac_list.txt ~/transcriptomes/pipeline_sample_data/sample_loci_shortened.txt ~/transcriptomes/pipeline_sample_data/sample_blastdb 4 Parallel 2
tfastq_assembly_master.py 7June2017 none ~/transcriptomes/pipeline_sample_data/Basellaceae/ t1_ ~/transcriptomes/shell/Basellaceae/ Bas ~/transcriptomes/pipeline_sample_data/exon_dbs/ none ~/transcriptomes/pipeline_sample_data/Basellaceae/Bas_list.txt ~/transcriptomes/pipeline_sample_data/sample_loci_shortened.txt ~/transcriptomes/pipeline_sample_data/sample_blastdb 4 Parallel 2
tfastq_assembly_master.py 7June2017 none ~/transcriptomes/pipeline_sample_data/Cactaceae/ t1_ ~/transcriptomes/shell/Cactaceae/ Cact ~/transcriptomes/pipeline_sample_data/exon_dbs/ none ~/transcriptomes/pipeline_sample_data/Cactaceae/Cact_list.txt ~/transcriptomes/pipeline_sample_data/sample_loci_shortened.txt ~/transcriptomes/pipeline_sample_data/sample_blastdb 4 Parallel 2

The following are what you will need to change to run it on your computer:
7June2017 is the date.
~/transcriptomes/pipeline_sample_data/ is the folder in which all of the sample data are found
~/transcriptomes/pipeline_sample_data/Anacampserotaceae/ ~/transcriptomes/pipeline_sample_data/Basellaceae/ and ~/transcriptomes/pipeline_sample_data/Cactaceae/ are the folders in which the fastq files are found.
~/transcriptomes/shell/Anacampserotaceae/ ~/transcriptomes/shell/Basellaceae/ and ~/transcriptomes/shell/Cactaceae/ are the scratch folders in which the analyses are actually run.  **Note that you must make this folder; it does not try to make it itself.**

The script to run (for me; in the case of Anacampserotaceae) is ~/transcriptomes/shell/Anacampserotaceae/7June2017_fastq_assembly_master.sh
Before running the script, you will need to make it executable (chmod u+x ~/transcriptomes/shell/Anacampserotaceae/7June2017_fastq_assembly_master.sh)


If it finishes correctly, there will be both sc_*.fa and sb3_*.out files in the ~/transcriptomes/shell/Anacampserotaceae/s2_final_Anac/spades_contigs/ folder.

**Part Two**

tcontig_classif_master.py 7June2017 none ~/transcriptomes/shell/ s2_final_ spades_ s ~/transcriptomes/shell/output/ shell ~/transcriptomes/pipeline_sample_data/sample_locus_list.txt ~/transcriptomes/pipeline_sample_data/sample_loci_shortened.txt ~/transcriptomes/pipeline_sample_data/sample_pardict.txt ~/transcriptomes/pipeline_sample_data/sample_outgroup_list.txt ~/transcriptomes/pipeline_sample_data/sample_group_list.txt none ~/transcriptomes/pipeline_sample_data/sample_ind_list_groups.txt ~/transcriptomes/pipeline_sample_data/backbone_trees/ new_ _al 8 6 Parallel

The following are what you will need to change to run it on your computer:
7June2017 is the date.
~/transcriptomes/pipeline_sample_data/ is the folder in which all of the sample data are found (including the fastq files).
~/transcriptomes/shell/ is the folder in which the folders with the sequences (e.g., Anacampserotaceae) are found.
~/transcriptomes/shell/output/ is the folder where the output should go; you must also make this folder first.
shell is the prefix for the output.  The fact that it is called shell is purely coincidental.  The thing that makes this script run in the shell is the Parallel at the end.

The script to run (for me) is /gpfs/scratch/ajm3/slurm/output/slurm_contig_classification_script1.sh

If it finishes correctly, there will be a RAxML_bestTree.shellfi_o150_* files for all of your gene families in the ~/transcriptomes/shell/output/shell_final/ folder.

**Part Three**

tgenefam_to_spptree_master.py 7June2017 none ~/transcriptomes/shell/output/ shell _allover150_al.fa _o150 ~/transcriptomes/pipeline_sample_data/sample_locus_list.txt ~/transcriptomes/pipeline_sample_data/sample_ind_list_groups_og.txt 4 ~/transcriptomes/pipeline_sample_data/RAxML_bestTree.sample_combined_22inds_4seqs Oryza_sativa_cds Parallel New

The following are what you will need to change to run it on your computer:
7June2017 is the date.
~/transcriptomes/pipeline_sample_data/ is the folder in which all of the sample data are found (including the fastq files).
~/transcriptomes/shell/output/ is the folder where the output should go.  It will make this folder.
shell is the prefix for the output, still purely coincidentally.
Part Three expects to have the Notung-2.8.1.6-beta.jar script in the scripts folder as well.  The remaining programs are presumed to be on the server and accessible normally.  If you get an error message associated with any of the programs, it is likely not present or called something else on your cluster/computer.

The script to run (for me) is ~/transcriptomes/shell/output/shell_Species_Tree_Making.sh
You will likely need to change the version number for the Notung script and the path to that script, for part three to work (in tcombpars_to_trees.py for round 1 and in tgenefam_to_spptree_master.py for round 2--search for Notung in either case and uncomment and change the correct line).

If it finishes correctly, there will be a folder called /home/abby/transcriptomes/shell/output/shell_prunedgenetrees2/shellpgts2_raxmlbs/ where you will be able to run astral (see below)
and a RAxML_bipartitions.shellpgtc2_combined_21inds_9seqs tree.

To run astral from within the /home/abby/transcriptomes/shell/output/shell_prunedgenetrees2/shellpgts2_raxmlbs/ folder:
java -jar ~/source/ASTRAL-master/Astral/astral.4.10.10.jar -i shellpgts2_bestTrees.tre -b bootstrap_filelist.txt -r 100 -o shellpgts2_astral_100bs.tre