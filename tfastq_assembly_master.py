#! /usr/bin/env python

#tfastq_assembly_master.py version 1.0 11 Jan. 2016 Abby Moore
#tfastq_assembly_master.py runs all of the scripts to get from fastq files to contigs that are classified by locus for tcontig_classif_master.py to deal with.

#depends on the following scripts:
'''
trans_fastq_to2_blast.py
tbaits_blastn_parse.py
tblast_to_fastq.py
tassembly_to_blast.py
tbaits_blastn_parse.py
tassembly_to_loci.py
(not done: cleanup script)

also:
spades
blast
'''

import sys

'''
tcontig_classif_master.py Date ScriptFolder SeqFolder DataFolder DataFolder1Pre DataFolder2Pre LLFileName ParDictFileName OGFileName GroupListFileName PolyListFileName IndListFileName AlFolder AlFilePre AlFilePost NCores NRounds
tcontig_classif_master.py 11Nov2015 none ~/transcriptomes/ none spades_amk ~/transcriptomes/sandbox/amk/ amk ~/transcriptomes/sandbox/amk/LocusList.txt ~/transcriptomes/general/pardict6.txt ~/transcriptomes/general/outgroup_list_new.txt ~/transcriptomes/general/Group_List_sand2.txt ~/transcriptomes/general/PolyList.txt ~/transcriptomes/general/Ln1_inds.txt ~/transcriptomes/general/phot_optimal_trees/ new_ _al 5 3
/gpfs/scratch/ajm3/eedwards/scripts/tcontig_classif_master.py 14Nov2015 /gpfs/scratch/ajm3/eedwards/scripts/ /gpfs/scratch/ajm3/eedwards/ s2_ spades_Ln1c2 s /gpfs/scratch/ajm3/eedwards/Ln1_comb2/ Ln1c2 /gpfs/scratch/ajm3/eedwards/general/Locus_List_all.txt /gpfs/scratch/ajm3/eedwards/general/loci_shortened_18S.txt /gpfs/scratch/ajm3/eedwards/general/pardict6.txt /gpfs/scratch/ajm3/eedwards/general/outgroup_list_new.txt /gpfs/scratch/ajm3/eedwards/general/Group_List_CB.txt /gpfs/scratch/ajm3/eedwards/general/PolyList.txt /gpfs/scratch/ajm3/eedwards/general/Ln1_inds.txt /gpfs/scratch/ajm3/eedwards/general/phot_optimal_trees/ new_ _al 4 3
'''

Usage = '''
tcontig_classif_master.py runs all of the script associated with contig
classification.
[the date to appear in the log file name]
[the folder in which the scripts are found, or none, if none]
[the folder in which the sequences folders are found]
[the main folder in which the data files are to be put, or "same", if it is the 
same as the sequence folder]
[the first prefix for the data subfolders, or none, if none]
[the second prefix for the data subfolders (after the group name), or none, if none]
[the file name for the list of loci]
[the file name with paralog [tab] locus]
[the file name for the pardict file]
[the file name for the outgroup information file]
[the file name for the list of groups (group[tab]group abbreviation)]
[the file name for the list of groups with polyploids]
[the file name for the file telling which individuals belong to which groups]
[the folder for the backbone alignments and trees]
[the prefix for the backbone alignments and trees]
[the suffix for the backbone alignments]
[number of cores for Oscar]
[number of rounds the analysis should be run]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 21:
	sys.exit("ERROR!  This script requires 20 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
Date = sys.argv[1]
ScriptFolder = sys.argv[2]
if ScriptFolder[-1] != "/":
	ScriptFolder += "/"
if ScriptFolder == "none/":
	ScriptFolder = ""
	ScriptFolder2 = "none"
else:
	ScriptFolder2 = ScriptFolder
SeqFolder = sys.argv[3]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
DataFolder = sys.argv[4]
if DataFolder == "same":
	DataFolder = SeqFolder
if DataFolder[-1] != "/":
	DataFolder += "/"
DataFolder1Pre = sys.argv[5]
DataFolder2Pre = sys.argv[6]
BSFilePre = sys.argv[7]
if BSFilePre == "none":
	BSFilePre = ""
OutFolder = sys.argv[8]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[9]
if OutFilePre == "none":
	OutFilePre = ""
LLFileName = sys.argv[10]
LPFileName = sys.argv[11]
ParDictFileName = sys.argv[12]
GroupListFileName = sys.argv[14]
IndListFileName = sys.argv[16]
BlastDBFileName = sys.argv[17]
NCores = sys.argv[20]

##gunzip the .fastq.gz files
##trans_fastq_to_2blast.py
#trans_fastq_to_2blast.py InFolder InFilePre BCFileName OutFolder OutFilePre BlastDB BlastFilePre
##blast script--in parallel
##tbaits_blastn_parse.py
#tbaits_blastn_parse.py InFolderSeq InFolderBl BCFileName LocusFileName AlFilePre BLFilePre OutFolder OutFilePre
##tblast_to_fastq.py
#tblast_to_fastq.py InFileName SeqFolder SeqFilePre OutFolder OutFilePre OutMode[together separate]
##spades--in parallel
##tassembly_to_blast.py
#tassembly_to_blast.py InFolder LocusFileName SeqFolder SeqFilePre BCFileName BlastDBName OutFolder OutFilePre
##blast--in parallel
##tbaits_blastn_parse.py
#tbaits_blastn_parse.py InFolderSeq InFolderBl BCFileName LocusFileName AlFilePre BLFilePre OutFolder OutFilePre
##spades--in parallel
##tassembly_to_loci.py
#tassembly_to_loci.py InFolder SeqFilePre LocusFileName OutFolder OutFilePre BlastFolder BlastDBPost BlastFilePre AMode [spades, masurca, minimo]
##blast--in parallel
##cleanup--write the script here, but maybe don't run it automatically?

#NoneDictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value, but if the value is "none",
#then the value is ""
#modified from DictFromFile from tbaits_intron_removal.py
def NoneDictFromFile(FileName):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]] = Line[1]
		if Line[1] == "none":
			TempDict[Line[0]] = ""
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is GroupDict

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
def CaptureColumn(FileName, ColNum):
	TempList = [ ]
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r').split('\t')
		TempList.append(Line[ColNum])
	InFile.close()
	print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is PolyList

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

GroupDict = NoneDictFromFile(GroupListFileName)
PolyList = CaptureColumn(PolyListFileName, 0)

#Letting the final scripts know how many rounds this has been run for.
RoundList = [ ]
for RoundNum in range(1,NRounds+1):
	RoundList.append(OutFilePre+"_round"+str(RoundNum)+"\t"+OutFilePre+"rnd"+str(RoundNum)+"\n")
OutFileName = OutFolder+OutFilePre+"_InFileList.txt"
OutFileWriting(OutFileName, RoundList)

OutScript = [ ]
Line = "#! /bin/bash\n#SBATCH -J "+OutFilePre+"\n#SBATCH -t 12:00:00\n#SBATCH -n "+NCores+"\n#SBATCH --mem="+str(int(NCores)*8)+"G\n"
Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\nmodule load mafft\nmodule load raxml\n"
OutScript.append(Line)
for Group in GroupDict:
	#this bit is very simple, so it doesn't need to be run in parallel.
	DataPath1 = DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre
	if (DataFolder1Pre == "none") and (GroupDict[Group] == ""):
		DataPath1 = DataFolder+Group+"/"+DataFolder2Pre
	DataPath = DataPath1+OutFilePre
	#preparing the folders and running tbaits_intron_removal.py
	Line = "mkdir "+DataPath+"_exons/\n"
	Line += "mkdir "+DataPath+"_paralogs/\n"
	Line += "mkdir "+DataPath+"_sequences/\n"
	OutScript.append(Line)
	#Line = ScriptFolder+"tbaits_intron_removal.py "+LPFileName+" "+BSFilePre+"b3_ "+BSFilePre+"c_ "+DataPath1+"contigs/ same "+DataPath+"_exons/ "+OutFilePre+"e1_ "+OGFileName+" "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+ScriptFolder2+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	#Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	#OutScript.append(Line)
for RoundNum in range(1,NRounds+1):
	GroupFileList = [ ]
	GroupFileListName = OutFolder+OutFilePre+"GroupFileList"+str(RoundNum)+".sh"
	for Group in GroupDict:
		OutScriptGroup = [ ]
		#setting up the file names
		OutScriptGroupName = OutFolder+OutFilePre+Group+"Script"+str(RoundNum)+".sh"
		GroupFileList.append("chmod u+x "+OutScriptGroupName+" && "+OutScriptGroupName+"\n")
		DataPath1 = DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre
		if (DataFolder1Pre == "none") and (GroupDict[Group] == ""):
			DataPath1 = DataFolder+Group+"/"+DataFolder2Pre
		DataPath = DataPath1+OutFilePre
		#two options for running the exon-analysis scripts:
		#the option for if the original exons have all been classified:
		if RoundNum != 1:
			Line = "chmod u+x "+DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh\n"
			Line += DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh\n"
			Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
			OutScriptGroup.append(Line)
		#and the option for if they haven't been
		'''
		Line = "chmod u+x "+DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
		OutScriptGroup.append(Line)
		'''
		if RoundNum != NRounds:
			#******tseq_placer_dup.py (merge)
			Line = ScriptFolder+"tseq_placer_dup.py "+LLFileName+" "+ParDictFileName+" "+DataPath+"_exons/ "+OutFilePre+"e"+str(RoundNum)+"_ none same "+OutFilePre+"e"+str(RoundNum)+"_ "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" 0.10 merge >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
		else:
			#********tseq_placer_dup.py (separate)
			Line = ScriptFolder+"tseq_placer_dup.py "+LLFileName+" "+OutFolder+OutFilePre+"_round"+str(RoundNum-1)+"/"+OutFilePre+"rnd"+str(RoundNum-1)+"_pardict_new.txt "+DataPath+"_exons/ "+OutFilePre+"e"+str(RoundNum)+"_ none same "+OutFilePre+"e"+str(RoundNum)+"_ "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+OutFolder+OutFilePre+"_round"+str(RoundNum-1)+"/ "+OutFilePre+"rnd"+str(RoundNum-1)+"_ _allseqs_al 0.10 separate >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
		Line += "chmod u+x "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Analysis_Script.sh\n"
		OutScriptGroup.append(Line)
		if RoundNum != NRounds:
			#******tcontigs_to_fixed_paralogs.py (intermediate)
			Line = ScriptFolder+"tcontigs_to_fixed_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+DataPath+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ "+DataPath+"_exons/ "+OutFilePre+"e"+str(RoundNum+1)+"_ "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ScriptFolder2+" "+OGFileName+" intermediate none 0 "+OutFilePre+"e"+str(RoundNum)+"_ >> "+OutFolder+Date+"_"+OutFilePre+".log\n\n"
		else:
			#********tcontigs_to_fixed_paralogs.py (final/finalpoly)
			if Group in PolyList:
				Line = ScriptFolder+"tcontigs_to_fixed_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+DataPath+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ same none "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ScriptFolder2+" "+OGFileName+" finalpoly "+DataFolder+Group+"/"+GroupDict[Group]+"_poly.txt 2.5 "+OutFilePre+"e"+str(RoundNum)+"_ >> "+OutFolder+Date+"_"+OutFilePre+".log\n\n"
			else:
				Line = ScriptFolder+"tcontigs_to_fixed_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+DataPath+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ same none "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ScriptFolder2+" "+OGFileName+" final none 0 "+OutFilePre+"e"+str(RoundNum)+"_ >> "+OutFolder+Date+"_"+OutFilePre+".log\n\n"
		OutScriptGroup.append(Line)
		OutFileWriting(OutScriptGroupName, OutScriptGroup)
	OutFileWriting(GroupFileListName, GroupFileList)
	Line = "cat "+GroupFileListName+" | parallel --jobs "+NCores+" --progress\n"
	OutScript.append(Line)
	##Combining everything:
	#*****tparalog_combiner.py
	#first, making the directory
	Line = "mkdir "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"\n"
	OutScript.append(Line)
	Line = "\n\n\n"+ScriptFolder+"tparalog_combiner.py "+GroupListFileName+" "+OGFileName+" "+DataFolder+" "+DataFolder1Pre+" "+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+ScriptFolder2+" "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ParDictFileName+" "+NCores+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	OutScript.append(Line)
	if RoundNum != NRounds:
		Line = "chmod u+x "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/"+OutFilePre+"rnd"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += OutFolder+OutFilePre+"_round"+str(RoundNum)+"/"+OutFilePre+"rnd"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += "chmod u+x "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/"+OutFilePre+"rnd"+str(RoundNum)+"_file_moving_script.sh\n"
		Line += OutFolder+OutFilePre+"_round"+str(RoundNum)+"/"+OutFilePre+"rnd"+str(RoundNum)+"_file_moving_script.sh\n\n\n"
		OutScript.append(Line)
	#*****tbaits_cleanup.py
	#Need to add this.
#Combining everything after the last round:
#*******tundivcontigs_combiner.py
#first, making the directory
Line = "mkdir "+OutFolder+OutFilePre+"_contigsplit\n"
OutScript.append(Line)
Line = ScriptFolder+"tundivcontigs_combiner.py "+GroupListFileName+" "+DataFolder+" "+DataFolder1Pre+" "+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"s"+str(NRounds)+"_ "+OutFolder+OutFilePre+"_round"+str(NRounds-1)+"/ "+OutFilePre+"rnd"+str(NRounds-1)+"_ _allseqs_al "+OutFolder+OutFilePre+"_contigsplit/ "+OutFilePre+"cs_ "+NCores+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
Line += "chmod u+x "+OutFolder+OutFilePre+"_contigsplit/"+OutFilePre+"cs_Redo_Contigs_Alignment.sh\n"
Line += OutFolder+OutFilePre+"_contigsplit/"+OutFilePre+"cs_Redo_Contigs_Alignment.sh\n"
OutScript.append(Line)
#*******tcontig_selection.py
Line = ScriptFolder+"tcontig_selection.py "+OutFolder+OutFilePre+"_contigsplit/"+OutFilePre+"cs_Contigs_to_Redo.txt "+GroupListFileName+" "+OutFolder+OutFilePre+"_contigsplit/ "+OutFilePre+"cs_ same "+OutFilePre+"cs_ "+DataFolder+" "+DataFolder1Pre+" "+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"s"+str(NRounds)+"_ "+IndListFileName+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
OutScript.append(Line)
#*******tparcomb_combiner.py
#first, making the directory
Line = "mkdir "+OutFolder+OutFilePre+"_combined\n"
OutScript.append(Line)
Line = ScriptFolder+"tparcomb_combiner.py "+OutFolder+" "+OutFolder+OutFilePre+"_InFileList.txt "+IndListFileName+" "+OGFileName+" "+OutFolder+OutFilePre+"_contigsplit/ "+OutFilePre+"cs_ "+OutFolder+OutFilePre+"_combined/ "+OutFilePre+"cb_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+ScriptFolder2+" "+ParDictFileName+" "+NCores+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
Line += "chmod u+x "+OutFolder+OutFilePre+"_combined/"+OutFilePre+"cb_Locus_Analysis_Script.sh\n"
Line += OutFolder+OutFilePre+"_combined/"+OutFilePre+"cb_Locus_Analysis_Script.sh\n"
Line += "\n\n"
OutScript.append(Line)
#*******tparcomb_final.py
Line = "mkdir "+OutFolder+OutFilePre+"_final\n"
Line += ScriptFolder+"tparcomb_final.py "+OutFolder+OutFilePre+"_combined/ "+OutFilePre+"cb_ "+OutFolder+OutFilePre+"_final/ "+OutFilePre+"fi_ "+ScriptFolder2+" "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+OGFileName+" "+IndListFileName+" "+NCores+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
OutScript.append(Line)
Line = "chmod u+x "+OutFolder+OutFilePre+"_final/"+OutFilePre+"fi_Locus_Analysis_Script.sh\n"
Line += OutFolder+OutFilePre+"_final/"+OutFilePre+"fi_Locus_Analysis_Script.sh\n"
OutScript.append(Line)

OutFileName = OutFolder+OutFilePre+"_contig_classification_script.sh"
OutFileWriting(OutFileName, OutScript)

#making the script to reclassify the ambiguous paralogs, if necessary

OutScript2 = [ ]
Line = "#! /bin/bash\n#SBATCH -J "+OutFilePre+"\n#SBATCH -t 4:00:00\n#SBATCH -n "+NCores+"\n#SBATCH --mem="+str(int(NCores)*4)+"G\n"
Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\nmodule load mafft\n"
OutScript2.append(Line)
Line = ScriptFolder+"tambig_seq_renaming.py "+OutFolder+OutFilePre+"_final/ "+OutFilePre+"fi_ "+OutFolder+OutFilePre+"_contigsplit/ "+OutFilePre+"cs_ "+OutFolder+OutFilePre+"_final/ "+OutFilePre+"fi_ "+ScriptFolder2+" "+OGFileName+" "+IndListFileName+" "+NCores+" "+OutFolder+OutFilePre+"_combined/"+OutFilePre+"cb_Ambig_Paralog_List.txt >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
Line += "chmod u+x "+OutFolder+OutFilePre+"_final/"+OutFilePre+"fi_ "+"seq_renaming_script1.sh\n"
Line += OutFolder+OutFilePre+"_final/"+OutFilePre+"fi_ "+"seq_renaming_script1.sh\n"
OutScript2.append(Line)
OutFileName2 = OutFolder+OutFilePre+"_Ambiguous_Sequence_Renaming.sh"
OutFileWriting(OutFileName2, OutScript2)
print("Once %s has run, look at the trees to determine whether any 'Ambig' sequences can really be classified into one of the paralogs." % (OutFileName))
print("If they can be, modify the file %s by adding [tab] actual paralog to each line, save it with the same file name, and run the script %s.\n" % (OutFolder+OutFilePre+"_combined/"+OutFilePre+"cb_Ambig_Paralog_List.txt", OutFileName2))
sys.stderr.write("Once %s has run, look at the trees to determine whether any 'Ambig' sequences can really be classified into one of the paralogs." % (OutFileName))
sys.stderr.write("If they can be, modify the file %s by adding [tab] actual paralog to each line, save it with the same file name, and run the script %s.\n" % (OutFolder+OutFilePre+"_combined/"+OutFilePre+"cb_Ambig_Paralog_List.txt", OutFileName2))

