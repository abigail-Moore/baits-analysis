#! /usr/bin/env python

#torig_spp_tree_master.py version 1.0 14 April 2016 Abby Moore
#This script writes a master script to make the original species tree by getting
#ndhF, rbcL, matK, and ITS sequences from the baits.  It takes the group list
#as input, so it makes separate scripts for each group.

import subprocess
import sys

'''
torig_spp_tree_master.py GroupDictFN InFolder SeqFolder ScriptPath LLFileName Mode
/users/ajm3/data/ajm3/scripts/torig_spp_tree_master.py 14April2016 /users/ajm3/data/ajm3/general/Group_List_Ln3_4_5.txt /gpfs/scratch/ajm3/eedwards/ /users/ajm3/data/ajm3/ /users/ajm3/data/ajm3/scripts/ /users/ajm3/data/ajm3/general/Locus_List_spptree.txt /users/ajm3/data/ajm3/general/spp_tree_seq_database Array
/users/ajm3/data/ajm3/scripts/torig_spp_tree_master.py 14April2016 /users/ajm3/data/ajm3/general/GroupList_Ln12.txt /gpfs/scratch/ajm3/eedwards/ /users/ajm3/data/ajm3/ /users/ajm3/data/ajm3/scripts/ /users/ajm3/data/ajm3/general/Locus_List_spptree.txt /users/ajm3/data/ajm3/general/spp_tree_seq_database Array
'''

Usage = '''
torig_spp_tree_master.py
This script makes a first guess species tree by finding ndhF, rbcL, matK, and 
ITS sequences from the baits, assembling them with spades, and taking the 
longest contig as the most probable correct one.
[date]
[file name for the group file: GroupName [tab] GroupAbbreviation]
[working folder (containing fasta sequence files)]
[folder containing fastq sequences]
[folder in which the scripts are found, or 'none', if it they are on the default
path]
[file with the list of loci]
[name of the blast database]
[folder with blast databases for the individual loci]
[mode, either Parallel or Array]
'''

ModeList = ['Parallel', 'Array']

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 10:
	sys.exit("ERROR!!  This script requires 9 additional arguments and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
Date = sys.argv[1]
GroupDictFN = sys.argv[2]
InFolder = sys.argv[3]
if InFolder[-1] != "/":
	InFolder += "/"
SeqFolder = sys.argv[4]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
ScriptPath = sys.argv[5]
if ScriptPath == "none":
	ScriptPath = ""
LLFileName = sys.argv[6]
BlastDBName = sys.argv[7]
BlastDBFolder = sys.argv[8]
if BlastDBFolder[-1] != "/":
	BlastDBFolder += "/"
Mode = sys.argv[9]
if Mode not in ModeList:
	sys.exit("ERROR!!  You wanted mode %s, but the mode must be one of the following: %s.\n" % (Mode, ", ".join(ModeList)))
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

#DictFromFile makes a dictionary from a tab-delimited file, where the keys are in columns KC, and
#the values are in column VC
#from tcontig_selection.py, which was a modified version of the function in from tbaits_intron_removal.py
def DictFromFile(FileName, KC, VC):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[KC]] = Line[VC]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is OutIndDict

#################################################################################

GroupDict = DictFromFile(GroupDictFN, 0, 1)#GroupDict[GroupName]=GroupPre

NumHrs = 1
#make a different script for each group
if Mode == "Array":
	SBatchList = [ ]
for Group in GroupDict:
	InFolderG = InFolder+Group+"/"
	SeqFolderG = SeqFolder+Group+"/"
	OutScript = ['#! /bin/bash\n#SBATCH -J spptree_'+GroupDict[Group]+'\n#SBATCH -t '+str(NumHrs)+':00:00\n#SBATCH -n 1\nmodule load blast\nmodule load spades\n\n']
	#make the new folder
	Line = 'mkdir '+InFolderG+'stb_'+GroupDict[Group]+'\n'
	#OutScript.append(Line)
	#run torig_spp_tree_blasting.py
	#two alternatives for this line
	#alternative for sequences that are broken up into smaller fasta files
	Line = ScriptPath+"torig_spp_tree_blasting.py "+InFolderG+" a1_ "+InFolderG+"a1_Num_Ind_Files_List.txt "+BlastDBName+" "+InFolderG+"stb_"+GroupDict[Group]+" stb_ Parallel 4 >> "+InFolderG+Date+"_spptree.log\n"
	Line += "chmod u+x "+InFolderG+'stb_'+GroupDict[Group]+'/stb_BlastScript1.sh\n'
	Line += InFolderG+'stb_'+GroupDict[Group]+'/stb_BlastScript1.sh\n'
	#OutScript.append(Line)
	#alternative for sequences that are not broken up into smaller fasta files (directly from gzipped fastq files)
	Line = "ls "+SeqFolderG+"*.gz | parallel gunzip {}\n"
	Line += ScriptPath+"trans_fastq_to_2blast.py "+SeqFolderG+" t1_ "+SeqFolderG+Group+"_ind_list.txt "+InFolderG+" a1_ "+BlastDBName+" stb_"+GroupDict[Group]+"/stb_ Parallel 4 None >> "+InFolderG+Date+"_spptree.log\n"
	Line += "chmod u+x "+InFolderG+'stb_'+GroupDict[Group]+'/stb_blast_script1.sh\n'
	Line += InFolderG+'stb_'+GroupDict[Group]+'/stb_blast_script1.sh\n'
	#OutScript.append(Line)
	#parse the blast output and run a spades assembly
	Line = ScriptPath+"tbaits_blastn_parse.py "+InFolderG+"stb_"+GroupDict[Group]+"/ "+InFolderG+"a1_Num_Ind_Files_List.txt "+LLFileName+" stb_ "+InFolderG+"stb_"+GroupDict[Group]+"/ stb_ >> "+InFolderG+Date+"_spptree.log\n"
	#OutScript.append(Line)
	#**********tblast_to_fastq.py
	Line = "mkdir "+InFolderG+"sts_"+GroupDict[Group]+"\n"
	Line += ScriptPath+"tblast_to_fastq.py "+InFolderG+"stb_"+GroupDict[Group]+"/stb_Seqs_to_Loci.txt "+SeqFolderG+" t1_ "+InFolderG+"sts_"+GroupDict[Group]+"/ sts_ separate >> "+InFolderG+Date+"_spptree.log\n"
	Line += "date >> "+InFolderG+Date+"_spptree.log\n"
	#OutScript.append(Line)
	#running spades a second time
	Line = "chmod u+x "+InFolderG+"sts_"+GroupDict[Group]+"/spades_script_separate1.sh\n"
	Line += InFolderG+"sts_"+GroupDict[Group]+"/spades_script_separate1.sh\n"
	Line += "date >> "+InFolderG+Date+"_spptree.log\n"
	#OutScript.append(Line)
	#***********tassembly_to_loci.py
	Line = "mkdir "+InFolderG+"sts_final_"+GroupDict[Group]+"\n"
	Line += "mkdir "+InFolderG+"sts_final_"+GroupDict[Group]+"/spades_contigs\n"
	Line += ScriptPath+"tassembly_to_loci.py "+InFolderG+"sts_"+GroupDict[Group]+"/ sts_ "+InFolderG+"sts_"+GroupDict[Group]+"/sts_Locus_List.txt "+InFolderG+"sts_final_"+GroupDict[Group]+"/spades_contigs/ stsc_ "+BlastDBFolder+" _spp_tree stsb_ spades >> "+InFolderG+Date+"_spptree.log\n"
	Line += "date >> "+InFolderG+Date+"_spptree.log\n"
	Line += "chmod u+x "+InFolderG+"sts_final_"+GroupDict[Group]+"/spades_contigs/stsb_BlastScript1.sh\n"
	Line += InFolderG+"sts_final_"+GroupDict[Group]+"/spades_contigs/stsb_BlastScript1.sh\n"
	#OutScript.append(Line)
	#********tbaits_intron_removal.py
	Line = ScriptPath+"tbaits_to_spptreeseqs.py "+LLFileName+" stsb_ stsc_ "+InFolderG+"sts_final_"+GroupDict[Group]+"/spades_contigs/ same same stsf_  >> "+InFolderG+Date+"_spptree.log\n"
	OutScript.append(Line)
	OutFileName = InFolderG+"st_analysis_master_"+GroupDict[Group]+".sh"
	OutFileWriting(OutFileName, OutScript)
	if Mode == "Array":
		OutLine = "sbatch "+OutFileName
		SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		SBatchList.append(SBatchOut.strip('\r').strip('\n').split(" ")[-1])
if Mode == "Array":
	print("%d separate jobs were submitted to analyze the groups.\n" % (len(SBatchList)))
	sys.stderr.write("%d separate jobs were submitted to analyze the groups.\n" % (len(SBatchList)))
elif Mode == "Parallel":
	print("%d separate scripts were written, with names such as %s, which must now be run separately.\n" % (len(GroupDict), OutFileName))
	sys.stderr.write("%d separate scripts were written, with names such as %s, which must now be run separately.\n" % (len(GroupDict), OutFileName))

