#! /usr/bin/env python

#ttree_tester.py version 1.0 24 June 2015 Abby Moore
#This script writes a bash script that analyzes all groups of spades results on
#the desired tree or group of trees and concatenates the results.

#ttree_tester.py

import sys

Usage = '''
ttree_tester.py is a script to analyze all of the spades output (or other types 
of contigs, once they have been turned into fasta files) to make alignments of
paralog sequences and concatenates the results.
ttree_tester.py
[file with list of loci]
[text file with the names of the different folders, with the name of the sequences
folder in the first column and the shortened name for the subfolders in the 2nd]
[tab-delimitted file with locus in the first column and outgroup sequence name 
in the second column]
[file with the new locus in the first column and the old locus/paralog in the 
second column]
[file with the locus name in the first column, the sequence name in the second 
column, and the paralog to which that sequence belongs in the third column]
[folder in which data files are found]
[type of alignment contigs: spades, minimo, or mazurca]
[folder in which scripts are found]
[folder in which backbone alignments are found]
[prefix of the backbone alignment files, or "none" if none]
[suffix of the backbone alignment files, or "none" if none]
[output folder name--not the full path; output folders will be made inside the sequence folders]
[prefix for output files, or "none" if none]
'''

ContigTypeList = ['spades', 'minimo', 'masurca']

#ttree_tester.py LLFileName GFileName OGFileName SLFileName PFileName DataFolder ContigType ScriptFolder AlFolder AlFilePre AlFilePost TreeFolder TreeFilePre OutFolder OutFilePre Date
#/users/ajm3/scratch/eedwards/ttree_tester.py /users/ajm3/scratch/eedwards/general/Locus_List_alaAT.txt /users/ajm3/scratch/eedwards/general/Group_List_ADLM.txt /users/ajm3/scratch/eedwards/general/outgroup_list_new.txt /users/ajm3/scratch/eedwards/general/loci_shortened.txt /users/ajm3/scratch/eedwards/general/pardict4.txt /users/ajm3/scratch/eedwards/ spades /users/ajm3/scratch/eedwards/ /users/ajm3/scratch/eedwards/general/phot_optimal_trees/ new_ _al /users/ajm3/scratch/eedwards/general/phot_optimal_trees/ new_ /users/ajm3/scratch/eedwards/combined_ADLM/alaAT/ alaAT_ 20150625

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 17:
	sys.exit("ERROR!  This script requires 16 additional arguments (sorry....) and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
LLFileName = sys.argv[1]
GFileName = sys.argv[2]
OGFileName = sys.argv[3]
SLFileName = sys.argv[4]
PFileName = sys.argv[5]
DataFolder = sys.argv[6]
if DataFolder[-1] != "/":
	DataFolder += "/"
ContigType = sys.argv[7]
if (ContigType in ContigTypeList) == False:
	sys.exit("ERROR! You requested the contig type %s, but the only available contig types are %s." % (ConitgType, ",".join(ContigTypeList)))
ScriptFolder = sys.argv[8]
if ScriptFolder[-1] != "/":
	ScriptFolder += "/"
AlFolder = sys.argv[9]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[10]
AlFilePost = sys.argv[11]
TreeFolder = sys.argv[12]
if TreeFolder[-1] != "/":
	TreeFolder += "/"
TreeFilePre = sys.argv[13]
OutFolder = sys.argv[14]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[15]
if OutFilePre == "none":
	OutFilePre = ""
if OutFilePre[-1] == "_":
	OutFilePre = OutFilePre[:-1]
Date = sys.argv[16]

#DictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
#from tbaits_intron_removal.py
def DictFromFile(FileName):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]] = Line[1]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict

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
	return

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
	#This is LocusList

if ContigType == "spades":
	DataFolder1Pre = "s2_"
	DataFolder2Pre = "spades_"
	DataFilePre = "s"
elif ContigType == "minimo":
	DataFolder1Pre = "min_"
	DataFolder2Pre = "minimo_"
	DataFilePre = "min"
elif ContigType == "masurca":
	DataFolder1Pre = "mas_"
	DataFolder2Pre = "masurca_"
	DataFilePre = "mas"

#getting the list of loci from the file
LocusList = CaptureColumn(LLFileName, 0)
#make the dictionary of the taxonomic groups and their shortened names
#GroupDict[Group] = GroupAbb
GroupDict = DictFromFile(GFileName)

#making the BlastFileList
#desired form: sb3_mdh.out	sc_mdh.fa
BlastFileList = [ ]
for Locus in LocusList:
	Line = DataFilePre+"b3_"+Locus+".out\t"+DataFilePre+"c_"+Locus+".fa\n"
	BlastFileList.append(Line)
OutFileName = OutFolder+OutFilePre+"_BlastFileList.txt"
OutFileWriting(OutFileName, BlastFileList)

OutScript = [ ]
Line = "#! /bin/bash\n#SBATCH -J "+OutFilePre+"\n#SBATCH -t 8:00:00\n#SBATCH -n 1\n#SBATCH --mem=16G\n"
Line += "mkdir "+OutFolder+"\n"
Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\nmodule load mafft\nmodule load raxml\n"
OutScript.append(Line)
for Group in GroupDict:
	DataPath =  DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre+OutFilePre
	DataPathOld = DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre
	Line = "mkdir "+DataPath+"_exons/\n"
	Line += ScriptFolder+"tbaits_intron_removal.py "+OutFolder+OutFilePre+"_BlastFileList.txt "+SLFileName+" "+DataFilePre+"c_ "+DataPathOld+"contigs/ same "+DataPath+"_exons/ "+OutFilePre+"e_ "+OGFileName+" "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+TreeFolder+" "+TreeFilePre+" "+ScriptFolder+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	Line += "chmod u+x "+DataPath+"_exons/"+OutFilePre+"e_Analysis_Script.sh\n"
	Line += DataPath+"_exons/"+OutFilePre+"e_Analysis_Script.sh\n"
	Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	OutScript.append(Line)
	Line = "mkdir "+DataPath+"_paralogs/\n"
	Line += ScriptFolder+"tseq_placer.py "+DataPath+"_exons/"+OutFilePre+"e_Locus_List.txt "+PFileName+" "+DataPath+"_exons/ "+OutFilePre+"e_ same "+OutFilePre+"e_ "+DataPath+"_paralogs/ "+OutFilePre+"p1_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	Line += "chmod u+x "+DataPath+"_paralogs/"+OutFilePre+"p1_Analysis_Script.sh\n"
	Line += DataPath+"_paralogs/"+OutFilePre+"p1_Analysis_Script.sh\n"
	Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	OutScript.append(Line)
	Line = "mkdir "+DataPath+"_sequences/\n"
	Line += ScriptFolder+"tcontigs_to_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p1_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p1_ "+DataPath+"_sequences/ "+OutFilePre+"s1_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+ScriptFolder+" >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	#But actually we want to use the other script instead of this one.
	#Line += "chmod u+x "+DataPath+"_sequences/"+OutFilePre+"s1_Analysis_Script.sh\n"
	#Line += DataPath+"_sequences/"+OutFilePre+"s1_Analysis_Script.sh\n"
	Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
	OutScript.append(Line)
	#************At this point there needs to be another program that analyzes the results of the previous stuff and makes some output including Locus_Paralog_List.txt!!!
Line = "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
Line += ScriptFolder+"tparalog_combiner.py "+GFileName+" "+OGFileName+" "+DataFolder+" "+DataFolder2Pre+" /"+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"1_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+OutFolder+" "+OutFilePre+"_ >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
Line += "chmod u+x "+OutFolder+OutFilePre+"Analysis_Script.sh\n"
Line += OutFolder+OutFilePre+"Analysis_Script.sh\n"
Line += "date >> "+OutFolder+Date+"_"+OutFilePre+".log\n"
OutScript.append(Line)
OutFileName = OutFolder+OutFilePre+"_Master_Analysis_Script.sh"
OutFileWriting(OutFileName, OutScript)
print("Now you need to make the script %s executable and execute it.\n" % (OutFileName))
sys.stderr.write("Now you need to make the script %s executable and execute it.\n" % (OutFileName))
