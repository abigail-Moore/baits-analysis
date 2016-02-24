#! /usr/bin/env python

#tfastq_assembly_master.py version 1.0 11 Jan. 2016 Abby Moore
#tfastq_assembly_master.py runs all of the scripts to get from fastq files to contigs that are classified by locus for tcontig_classif_master.py to deal with.

#depends on the following scripts:
'''
trans_fastq_to2_blast.py
tbaits_blastn_parse.py
tblast_to_fastq.py
tassembly_to_blast.py
tassembly_to_loci.py

also:
spades
blast
parallel
'''

import sys

'''
tfastq_assembly_master.py Date ScriptFolder SeqFolder SeqFilePre DataFolder GroupPre BlastDBFolder BlastDBPost IndListFileName LLFileName BlastDBFileName NCores
'''

Usage = '''
tfastq_assembly_master.py runs all of the scripts associated with read assembly
into contigs.
[the date to appear in the log file name]
[the folder in which the scripts are found, or none, if none]
[the folder in which the sequences are found]
[the prefix for the sequence files, or "none", if none]
[the main folder in which the data files are to be put, or "same", if it is the 
same as the sequence folder]
[the abbreviated name for this group of sequences]
[the folder with the blast databases for the individual gene families]
[the suffix for the blast database names, or "none", if none]
[file with list of individuals, one individual per line]
[file with list of loci, one locus per line]
[the main blast database file containing all gene families]
[number of cores for Oscar]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 13:
	sys.exit("ERROR!  This script requires 12 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
Date = sys.argv[1]
ScriptFolder = sys.argv[2]
if ScriptFolder[-1] != "/":
	ScriptFolder += "/"
if ScriptFolder == "none/":
	ScriptFolder = ""
SeqFolder = sys.argv[3]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
SeqFilePre = sys.argv[4]
if SeqFilePre == "none":
	SeqFilePre2 = ""
else:
	SeqFilePre2 = SeqFilePre
DataFolder = sys.argv[5]
if DataFolder == "same":
	DataFolder = SeqFolder
elif DataFolder[-1] != "/":
	DataFolder += "/"
GroupPre = sys.argv[6]
BlastDBFolder = sys.argv[7]
if BlastDBFolder[-1] != "/":
	BlastDBFolder += "/"
BlastDBPost = sys.argv[8]
IndListFileName = sys.argv[9]
LLFileName = sys.argv[10]
BlastDBFileName = sys.argv[11]
NCores = sys.argv[12]

################################################################################

#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
#from tbaits_intron_removal.py
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	print("Output file %s written.\n" % (FileName))
	sys.stderr.write("Output file %s written.\n" % (FileName))

##################################################################################

OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'fastq_assembly\n#SBATCH -t 48:00:00\n#SBATCH -n '+NCores+'\n#SBATCH --mem='+str(int(NCores)*8)+'G\n\n']
Line = "module load blastn\nmodule load spades\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
##gunzip the .fastq.gz files
Line = "cat "+IndListFileName+" | parallel gunzip "+SeqFolder+SeqFilePre2+"{}_R1.fastq.gz "+SeqFolder+SeqFilePre2+"{}_R2.fastq.gz\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
##################################ROUND ONE######################################
#*********trans_fastq_to_2blast.py
Line = "mkdir "+DataFolder+"\n"
Line += "mkdir "+DataFolder+"b1_"+GroupPre+"\n"
Line += ScriptFolder+"trans_fastq_to_2blast.py "+SeqFolder+" "+SeqFilePre+" "+IndListFileName+" "+DataFolder+" a1_ "+BlastDBFileName+" b1_"+GroupPre+"/b1_ "+NCores+" >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running blast
Line = "chmod u+x "+DataFolder+"b1_"+GroupPre+"/b1_blast_script1.sh\n"
Line += DataFolder+"b1_"+GroupPre+"/b1_blast_script1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#*************tbaits_blastn_parse.py
Line = ScriptFolder+"tbaits_blastn_parse.py "+DataFolder+" "+DataFolder+"b1_"+GroupPre+"/ "+IndListFileName+" "+LLFileName+" a1_ b1_ "+DataFolder+"b1_"+GroupPre+"/ b1_ >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#***********tblast_to_fastq.py
Line = "mkdir "+DataFolder+"s1_"+GroupPre+"\n"
Line += ScriptFolder+"tblast_to_fastq.py "+DataFolder+"b1_"+GroupPre+"/b1_Seqs_to_Loci.txt "+SeqFolder+" "+SeqFilePre+" "+DataFolder+"s1_"+GroupPre+"/ s1_ together >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running spades
Line = "chmod u+x "+DataFolder+"s1_"+GroupPre+"/spades_script_together1.sh\n"
Line += DataFolder+"s1_"+GroupPre+"/spades_script_together1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
##***********tassembly_to_blast.py
Line = "mkdir "+DataFolder+"b2_"+GroupPre+"\n"
Line += ScriptFolder+"tassembly_to_blast.py "+DataFolder+"s1_"+GroupPre+"/ "+DataFolder+"s1_"+GroupPre+"/s1_Locus_List.txt "+DataFolder+" a1_ "+IndListFileName+" "+BlastDBFileName+" "+DataFolder+"b2_"+GroupPre+"/ b2_ >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)

##################################ROUND TWO######################################
#running blast a second time
Line = "chmod u+x "+DataFolder+"b2_"+GroupPre+"/b2_BlastScript1.sh\n"
Line += DataFolder+"b2_"+GroupPre+"/b2_BlastScript1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#**********tbaits_blastn_parse.py
Line = ScriptFolder+"tbaits_blastn_parse.py "+DataFolder+" "+DataFolder+"b2_"+GroupPre+"/ "+IndListFileName+" "+LLFileName+" a1_ b2_ "+DataFolder+"b2_"+GroupPre+"/ b2_ >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#**********tblast_to_fastq.py
Line = "mkdir "+DataFolder+"s2_"+GroupPre+"\n"
Line += ScriptFolder+"tblast_to_fastq.py "+DataFolder+"b2_"+GroupPre+"/b2_Seqs_to_Loci.txt "+SeqFolder+" "+SeqFilePre+" "+DataFolder+"s2_"+GroupPre+"/ s2_ separate >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running spades a second time
Line = "chmod u+x "+DataFolder+"s2_"+GroupPre+"/spades_script_separate1.sh\n"
Line += DataFolder+"s2_"+GroupPre+"/spades_script_separate1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#***********tassembly_to_loci.py
Line = "mkdir "+DataFolder+"s2_"+GroupPre+"_spades_contigs\n"
Line += ScriptFolder+"tassembly_to_loci.py "+DataFolder+"s2_"+GroupPre+"/ s2_ "+DataFolder+"s2_"+GroupPre+"/s2_Locus_List.txt "+DataFolder+"s2_"+GroupPre+"_spades_contigs/ sc_ "+BlastDBFolder+" "+BlastDBPost+" sb3_ spades >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running blast a third time to find the exons
Line = "chmod u+x "+DataFolder+"s2_"+GroupPre+"_spades_contigs/sb3_BlastScript1.sh\n"
Line += DataFolder+"s2_"+GroupPre+"_spades_contigs/sb3_BlastScript1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)

############################CLEANUP###############################################
##cleanup--write the script here, but maybe don't run it automatically?
CleanupScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'assembly_cleanup\n#SBATCH -t 4:00:00\n#SBATCH -n 1\n\n']
#rezip original sequence files
Line = "cat "+IndListFileName+" | parallel gzip "+SeqFolder+SeqFilePre2+"{}_R1.fastq "+SeqFolder+SeqFilePre2+"{}_R2.gz\n"
#remove fasta files
Line += "rm "+DataFolder+"a1_*.fa\n"
CleanupScript.append(Line)
##first round:
#remove first blast folder
Line = "rm -r "+DataFolder+"b1_"+GroupPre+"/\n"
#remove first spades folder 
Line += "rm -r "+DataFolder+"s1_"+GroupPre+"/\n"
CleanupScript.append(Line)
##second round:
#remove second blast folder 
Line = "rm -r "+DataFolder+"b2_"+GroupPre+"/\n"
#remove second spades folder
Line += "rm -r "+DataFolder+"s2_"+GroupPre+"/\n"
CleanupScript.append(Line)


OutFileName = DataFolder+Date+"_fastq_assembly_master.sh"
OutFileWriting(OutFileName, OutScript)
OutFileName = DataFolder+Date+"_fastq_cleanup.sh"
OutFileWriting(OutFileName, CleanupScript)
