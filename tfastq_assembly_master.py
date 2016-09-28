#! /usr/bin/env python

#tfastq_assembly_master.py version 1.1 14 March 2016 Abby Moore
#tfastq_assembly_master.py runs all of the scripts to get from fastq files to contigs that are classified by locus for tcontig_classif_master.py to deal with.
#version 1.0 11 Jan. 2016
#version 1.1 14 March 2016 modified to make job arrays for blast instead of running in parallel

#depends on the following scripts:
'''
trans_fastq_to_2blast.py
tbaits_blastn_parse.py
tblast_to_fastq.py
tassembly_to_blast.py
tassembly_to_loci.py
tblast_to_fasta_qual.py

also:
spades
blast
parallel
minimo

in Array mode:
Slurm
'''

import sys

'''
tfastq_assembly_master.py Date ScriptFolder SeqFolder SeqFilePre DataFolder GroupPre BlastDBFolder BlastDBPost IndListFileName LLFileName BlastDBFileName NCores Mode[Parallel, Array]
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
[mode: either Parallel (for desktops or clusters) or Array (for Slurm)]
[blast rounds: 1 or 2, for cases in which there is a good database of sequences
already (1) or when the baits are the only database for this group (2)]
'''

ModeList = ['Parallel', 'Array']

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 15:
	sys.exit("ERROR!  This script requires 14 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
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
Mode = sys.argv[13]
if Mode not in ModeList:
	sys.exit("ERROR!!  You wanted mode %s, but the mode can only be one of the following: %s.\n %s" % (Mode, ", ".join(ModeList), Usage))
BlastRounds = int(sys.argv[14])
if BlastRounds not in [1, 2]:
	sys.exit("ERROR!!  You wanted %d blast rounds, but there can only be 1 or 2.\n%s" % (BlastRounds, Usage))

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

if Mode == "Parallel":
	NumHrs = 48
elif Mode == "Array":
	NumHrs = 2
OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_fastq_assembly\n#SBATCH -t '+str(NumHrs)+':00:00\n#SBATCH -n '+NCores+'\n\n\n']#SBATCH --mem='+str(int(NCores)*8)+'G\n\n']
Line = "module load blast\nmodule load spades\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
##gunzip the .fastq.gz files
Line = "cat "+IndListFileName+" | parallel gunzip "+SeqFolder+SeqFilePre2+"{}_R1.fastq.gz "+SeqFolder+SeqFilePre2+"{}_R2.fastq.gz\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
##################################ROUND ONE######################################
#*********trans_fastq_to_2blast.py
Line = "mkdir "+DataFolder+"\n"
if BlastRounds == 1:
	Line += "mkdir "+DataFolder+"b2_"+GroupPre+"\n"
	Line += ScriptFolder+"trans_fastq_to_2blast.py "+SeqFolder+" "+SeqFilePre+" "+IndListFileName+" "+DataFolder+" a1_ "+BlastDBFileName+" b2_"+GroupPre+"/b2_ "+Mode+" "+NCores+" "+DataFolder+"fastq_assembly3.sh >> "+DataFolder+Date+"_fastq_assembly.log\n"
	OutScript.append(Line)
	if Mode == "Array":
		OutFileName = DataFolder+"fastq_assembly1.sh"
		OutFileWriting(OutFileName, OutScript)
		print("Since you are running in array mode, the file you will need to execute is %s.  It should then execute the remaining files.\n" % (OutFileName))
		sys.stderr.write("Since you are running in array mode, the file you will need to execute is %s.  It should then execute the remaining files.\n" % (OutFileName))
		OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_fastq_assembly2\n#SBATCH -t '+str(NumHrs*2)+':00:00\n#SBATCH -n '+NCores+'\n\n\nmodule load blast\nmodule load spades\n']
	#running blast
	elif Mode == "Parallel":
		Line = "chmod u+x "+DataFolder+"b2_"+GroupPre+"/b2_blast_script1.sh\n"
		Line += DataFolder+"b2_"+GroupPre+"/b2_blast_script1.sh\n"
		OutScript.append(Line)
elif BlastRounds == 2:
	Line += "mkdir "+DataFolder+"b1_"+GroupPre+"\n"
	Line += ScriptFolder+"trans_fastq_to_2blast.py "+SeqFolder+" "+SeqFilePre+" "+IndListFileName+" "+DataFolder+" a1_ "+BlastDBFileName+" b1_"+GroupPre+"/b1_ "+Mode+" "+DataFolder+"fastq_assembly2.sh >> "+DataFolder+Date+"_fastq_assembly.log\n"
	OutScript.append(Line)
	if Mode == "Array":
		OutFileName = DataFolder+"fastq_assembly1.sh"
		OutFileWriting(OutFileName, OutScript)
		print("Since you are running in array mode, the file you will need to execute is %s.  It should then execute the remaining files.\n" % (OutFileName))
		sys.stderr.write("Since you are running in array mode, the file you will need to execute is %s.  It should then execute the remaining files.\n" % (OutFileName))
		OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_fastq_assembly2\n#SBATCH -t '+str(NumHrs*4)+':00:00\n#SBATCH -n '+NCores+'\n#SBATCH --mem='+str(int(NCores)*8)+'G\n\n\nmodule load blast\nmodule load spades\n']
	#running blast
	elif Mode == "Parallel":
		Line = "chmod u+x "+DataFolder+"b1_"+GroupPre+"/b1_blast_script1.sh\n"
		Line += DataFolder+"b1_"+GroupPre+"/b1_blast_script1.sh\n"
		OutScript.append(Line)
	#*************tbaits_blastn_parse.py
	Line = "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
	Line += ScriptFolder+"tbaits_blastn_parse.py "+DataFolder+"b1_"+GroupPre+"/ "+DataFolder+"a1_Num_Ind_Files_List.txt "+LLFileName+" b1_ "+DataFolder+"b1_"+GroupPre+"/ b1_ >> "+DataFolder+Date+"_fastq_assembly.log\n"
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
	Line += ScriptFolder+"tassembly_to_blast.py "+DataFolder+"s1_"+GroupPre+"/ "+DataFolder+"s1_"+GroupPre+"/s1_Locus_List.txt "+DataFolder+" a1_ "+DataFolder+"a1_Num_Ind_Files_List.txt "+BlastDBFileName+" "+DataFolder+"b2_"+GroupPre+"/ b2_ "+Mode+" "+DataFolder+"fastq_assembly3.sh >> "+DataFolder+Date+"_fastq_assembly.log\n"
	OutScript.append(Line)

##################################ROUND TWO######################################
	#running blast a second time
	if Mode == "Array":
		OutFileName = DataFolder+"fastq_assembly2.sh"
		OutFileWriting(OutFileName, OutScript)
		print("Since you are running in array mode, the file %s will be automatically executed.  It should then execute the remaining files.\n" % (OutFileName))
		sys.stderr.write("Since you are running in array mode, the file %s will be automatically executed.  It should then execute the remaining files.\n" % (OutFileName))
		OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_fastq_assembly3\n#SBATCH -t '+str(NumHrs*8)+':00:00\n#SBATCH -n '+NCores+'\n\n\nmodule load blast\nmodule load spades\n']
	elif Mode == "Parallel":
		Line = "chmod u+x "+DataFolder+"b2_"+GroupPre+"/b2_BlastScript1.sh\n"
		Line += DataFolder+"b2_"+GroupPre+"/b2_BlastScript1.sh\n"
		OutScript.append(Line)
#**********tbaits_blastn_parse.py
Line = "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += ScriptFolder+"tbaits_blastn_parse.py "+DataFolder+"b2_"+GroupPre+"/ "+DataFolder+"a1_Num_Ind_Files_List.txt "+LLFileName+" b2_ "+DataFolder+"b2_"+GroupPre+"/ b2_ >> "+DataFolder+Date+"_fastq_assembly.log\n"
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
Line = "mkdir "+DataFolder+"s2_final_"+GroupPre+"\n"
Line += "mkdir "+DataFolder+"s2_final_"+GroupPre+"/spades_contigs\n"
Line += ScriptFolder+"tassembly_to_loci.py "+DataFolder+"s2_"+GroupPre+"/ s2_ "+DataFolder+"s2_"+GroupPre+"/s2_Locus_List.txt "+DataFolder+"s2_final_"+GroupPre+"/spades_contigs/ sc_ "+BlastDBFolder+" "+BlastDBPost+" sb3_ spades >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running blast a third time to find the exons
Line = "chmod u+x "+DataFolder+"s2_final_"+GroupPre+"/spades_contigs/sb3_BlastScript1.sh\n"
Line += DataFolder+"s2_final_"+GroupPre+"/spades_contigs/sb3_BlastScript1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
if Mode == "Array":
	OutFileName = DataFolder+"fastq_assembly3.sh"
	OutFileWriting(OutFileName, OutScript)
	print("Since you are running in array mode, the file %s will also be automatically executed.\n" % (OutFileName))
	sys.stderr.write("Since you are running in array mode, the file %s will also be automatically executed.\n" % (OutFileName))
elif Mode == "Parallel":
	OutFileName = DataFolder+Date+"_fastq_assembly_master.sh"
	OutFileWriting(OutFileName, OutScript)

############################CLEANUP###############################################
##cleanup--write the script here, but maybe don't run it automatically?
CleanupScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_assembly_cleanup\n#SBATCH -t 4:00:00\n#SBATCH -n 1\n\n']
#remove fasta files
Line = "rm "+DataFolder+"a1_*.fa\n"
CleanupScript.append(Line)
if BlastRounds == 2:
	##first round:
	#remove first blast folder
	Line = "rm -r "+DataFolder+"b1_"+GroupPre+"/\n"
	#remove first spades folder 
	Line += "rm -r "+DataFolder+"s1_"+GroupPre+"/\n"
	CleanupScript.append(Line)
##second round:
#remove the blast output files from the second blast folder, while leaving the .txt files and the blast database 
Line = "rm "+DataFolder+"b2_"+GroupPre+"/*.out\n"
#remove second spades folder
Line += "rm -r "+DataFolder+"s2_"+GroupPre+"/\n"
CleanupScript.append(Line)
OutFileName = DataFolder+Date+"_fastq_cleanup.sh"
OutFileWriting(OutFileName, CleanupScript)

###########################MINIMO SCRIPT##########################################
OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_minimo\n#SBATCH -t 24:00:00\n#SBATCH -n '+NCores+'\n\n']
Line = "module load amos\nmodule load blast\n"
Line += "mkdir "+DataFolder+"min_"+GroupPre+"/\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#preparing the files for minimo
Line = ScriptFolder+"tblast_to_fasta_qual.py "+DataFolder+"b2_"+GroupPre+"/b2_Seqs_to_Loci.txt "+SeqFolder+" "+SeqFilePre+" "+DataFolder+"min_"+GroupPre+"/ min_ separate "+NCores+" >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running minimo
Line = "chmod u+x "+DataFolder+"min_"+GroupPre+"/minimo_script_separate1.sh\n"
Line += DataFolder+"min_"+GroupPre+"/minimo_script_separate1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#***********tassembly_to_loci.py
Line = "mkdir "+DataFolder+"min_final_"+GroupPre+"\n"
Line += "mkdir "+DataFolder+"min_final_"+GroupPre+"/minimo_contigs\n"
Line += ScriptFolder+"tassembly_to_loci.py "+DataFolder+"min_"+GroupPre+"/ min_ "+DataFolder+"min_"+GroupPre+"/min_Locus_List.txt "+DataFolder+"min_final_"+GroupPre+"/minimo_contigs/ minc_ "+BlastDBFolder+" "+BlastDBPost+" minb3_ minimo >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running blast a third time to find the exons
Line = "chmod u+x "+DataFolder+"min_final_"+GroupPre+"/minimo_contigs/minb3_BlastScript1.sh\n"
Line += DataFolder+"min_final_"+GroupPre+"/minimo_contigs/minb3_BlastScript1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
OutFileName = DataFolder+Date+"_minimo_assembly.sh"
OutFileWriting(OutFileName, OutScript)
print("The script for running minimo is %s, but it needs to be run separately from the main script.\n" % (OutFileName))
sys.stderr.write("The script for running minimo is %s, but it needs to be run separately from the main script.\n" % (OutFileName))

##########################SSAKE SCRIPT###########################################
OutScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_ssake\n#SBATCH -t 4:00:00\n#SBATCH -n '+NCores+'\n\n']
Line = "module load blast\n"
Line += "mkdir "+DataFolder+"ssa_"+GroupPre+"/\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#preparing the files for ssake
Line = ScriptFolder+"tblast_to_ssake.py "+DataFolder+"b2_"+GroupPre+"/b2_Seqs_to_Loci.txt "+SeqFolder+" "+SeqFilePre+" "+DataFolder+"ssa_"+GroupPre+"/ ssa_ separate "+NCores+" "+ScriptFolder+" >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running ssake
Line = "chmod u+x "+DataFolder+"ssa_"+GroupPre+"/ssake_script_separate1.sh\n"
Line += DataFolder+"ssa_"+GroupPre+"/ssake_script_separate1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#***********tassembly_to_loci.py
Line = "mkdir "+DataFolder+"ssa_final_"+GroupPre+"\n"
Line += "mkdir "+DataFolder+"ssa_final_"+GroupPre+"/ssake_contigs\n"
#************tassembly_to_loci.py still needs to be changed to incorporate ssake
Line += ScriptFolder+"tassembly_to_loci.py "+DataFolder+"ssa_"+GroupPre+"/ ssa_ "+DataFolder+"ssa_"+GroupPre+"/ssa_Locus_List.txt "+DataFolder+"ssa_final_"+GroupPre+"/ssake_contigs/ ssac_ "+BlastDBFolder+" "+BlastDBPost+" ssab3_ ssake >> "+DataFolder+Date+"_fastq_assembly.log\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
#running blast a third time to find the exons
Line = "chmod u+x "+DataFolder+"ssa_final_"+GroupPre+"/ssake_contigs/ssab3_BlastScript1.sh\n"
Line += DataFolder+"ssa_final_"+GroupPre+"/ssake_contigs/ssab3_BlastScript1.sh\n"
Line += "date >> "+DataFolder+Date+"_fastq_assembly.log\n"
OutScript.append(Line)
OutFileName = DataFolder+Date+"_ssake_assembly.sh"
OutFileWriting(OutFileName, OutScript)
print("The script for running ssake is %s, but it needs to be run separately from the main script.\n" % (OutFileName))
sys.stderr.write("The script for running ssake is %s, but it needs to be run separately from the main script.\n" % (OutFileName))

#########################Second cleanup script###################################
##cleanup--write the script here, but don't run it automatically
CleanupScript = ['#! /bin/bash\n#SBATCH -J '+GroupPre+'_min_ssake_cleanup\n#SBATCH -t 4:00:00\n#SBATCH -n 1\n\n']
#rezip original sequence files
Line = "cat "+IndListFileName+" | parallel gzip "+SeqFolder+SeqFilePre2+"{}_R1.fastq "+SeqFolder+SeqFilePre2+"{}_R2.gz\n"
CleanupScript.append(Line)
#remove minimo folder
Line += "rm -r "+DataFolder+"min_"+GroupPre+"/\n"
CleanupScript.append(Line)
#remove ssake folder
Line += "rm -r "+DataFolder+"ssa_"+GroupPre+"/\n"
CleanupScript.append(Line)
OutFileName = DataFolder+Date+"_min_ssake_cleanup.sh"
OutFileWriting(OutFileName, CleanupScript)
