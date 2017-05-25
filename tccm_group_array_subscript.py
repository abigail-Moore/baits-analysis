#! /usr/bin/env python

#tccm_group_array_subscript.py version 1.0 23 March 2016
#This script is a subscript to be called by tcontig_classif_master.py to
#write arrays to classify sequences of the individual groups

#depends on:
'''
tseq_placer_dup.py
tcontigs_to_fixed_paralogs.py
tparalog_combiner.py

mafft
raxml
fasta_to_phylip.py
'''

from collections import defaultdict
import sys
import subprocess #We want to be able to talk to the command line.

#tccm_group_array_subscript.py ScriptFolder DataFolder DataFolder1Pre DataFolder2Pre OutFolder OutFilePre LLFileName ParDictFileName OGFileName GroupListFileName PolyListFileName IndListFileName AlFolder AlFilePre AlFilePost NRounds RoundNum

Usage = '''
tccm_group_array_subscript.py is a script to be called by tcontig_classif_master.py 
to run the scripts for the individual groups as an array.
[the folder in which the scripts are found, or none, if none]
[the main folder in which the data files are found]
[the first prefix for the data subfolders, or none, if none]
[the second prefix for the data subfolders (after the group name), or none, if none]
[the output folder]
[the prefix for the output files and folders, or none, if none]
[the file name for the list of loci]
[the file name for the pardict file]
[the file name for the outgroup information file]
[the file name for the list of groups (group[tab]group abbreviation)]
[the file name for the list of groups with polyploids]
[the file name for the file telling which individuals belong to which groups]
[the folder for the backbone alignments and trees]
[the prefix for the backbone alignments and trees]
[the suffix for the backbone alignments]
[number of rounds the analysis should be run]
[the current round number]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 18:
	sys.exit("ERROR!  This script requires 18 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
ScriptFolder = sys.argv[1]
if ScriptFolder[-1] != "/":
	ScriptFolder += "/"
if ScriptFolder == "none/":
	ScriptFolder = ""
	ScriptFolder2 = "none"
else:
	ScriptFolder2 = ScriptFolder
DataFolder = sys.argv[2]
if DataFolder[-1] != "/":
	DataFolder += "/"
DataFolder1Pre = sys.argv[3]
DataFolder2Pre = sys.argv[4]
OutFolder = sys.argv[5]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[6]
if OutFilePre == "none":
	OutFilePre = ""
LLFileName = sys.argv[7]
ParDictFileName = sys.argv[8]
OGFileName = sys.argv[9]
GroupListFileName = sys.argv[10]
PolyListFileName = sys.argv[11]
IndListFileName = sys.argv[12]
AlFolder = sys.argv[13]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[14]
AlFilePost = sys.argv[15]
NRounds = int(sys.argv[16])
RoundNum = int(sys.argv[17])

#################################################################################

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

#ListDictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value.  Each key has multiple values, and they
#will be made into a list
#modified from tparalaog_combiner.py
def ListDictFromFile(FileName, KeyCol, ListCol):
	TempDict = defaultdict(list)
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[KeyCol]].append(Line[ListCol])
	InFile.close()
	for Key in TempDict:
		ListTemp = TempDict[Key]
		ListTemp = list(set(ListTemp))
		TempDict[Key] = ListTemp
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is GroupIndDict


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

#OrderedListMaking makes a list in which objects are arranged in descending 
#(if Rev == True) or ascending (if Rev == False) order
#It takes a dictionary with DictTemp[Thing] = Size
#from tblast_to_fasta_qual.py
def OrderedListMaking(DictTemp, Rev):
	DictbySize = defaultdict(list)
	for Thing in DictTemp:
		Size = DictTemp[Thing]
		DictbySize[Size].append(Thing)
	ListTemp = [ ]#list of Things in descending order of size
	for Size in sorted(DictbySize.keys(), reverse = Rev):
		ListTemp += DictbySize[Size]
	return(ListTemp)
	#This is GroupList


#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
#from tbaits_intron_removal.py
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	#print("Output file %s written.\n" % (FileName))
	#sys.stderr.write("Output file %s written.\n" % (FileName))


#################################################################################

#reading information from files:
GroupDict = NoneDictFromFile(GroupListFileName)
if PolyListFileName != "none":
	PolyList = CaptureColumn(PolyListFileName, 0)
else:
	PolyList = [ ]
GroupIndDict = ListDictFromFile(IndListFileName, 0, 1)
#sort the groups by the number of individuals
GroupNumDict = { }
for Group in GroupIndDict:
	GroupNumDict[Group] = len(GroupIndDict[Group])
GroupList = OrderedListMaking(GroupNumDict, True)


#making group files, starting with the groups with the most individuals

OutScriptGroupName = OutFolder+OutFilePre+Group+"_"+str(RoundNum)+".sh"

SBatchList = [ ]
for Group in GroupList:
	NumInds = GroupNumDict[Group]
	if (NumInds > 15):
		NCores = 4
	else:
		NCores = 2
	OutScriptGroup = [ "#! /bin/bash\n#SBATCH -J "+OutFilePre+Group+"_"+str(RoundNum)+"\n"]
	Line = "#SBATCH -t 12:00:00\n#SBATCH -n "+str(NCores)+"\n"
	Line += "#SBATCH --mem="+(str(NCores*8))+"GB\n"
	Line += "module load mafft\nmodule load raxml\n"
	Line += "date >> "+OutFolder+OutFilePre+Group+".log\n"
	OutScriptGroup.append(Line)
	#setting up the file names
	DataPath1 = DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre
	if (DataFolder1Pre == "none") and (GroupDict[Group] == ""):
		DataPath1 = DataFolder+Group+"/"+DataFolder2Pre
	DataPath = DataPath1+OutFilePre
	#two options for running the exon-analysis scripts:
	#the option for if the original exons have all been classified:
	'''
	if RoundNum != 1:
		Line = "chmod u+x "+DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh\n"
		Line += "date >> "+OutFolder+OutFilePre+Group+".log\n"
		OutScriptGroup.append(Line)
	'''
	#and the option for if they haven't been
	##***************This script needs to be rewritten so that it can be parallelized (meaning that all of the commands dealing with the same alignment need to be on the same line*********
	Line = "cat "+DataPath+"_exons/"+OutFilePre+"e"+str(RoundNum)+"_Analysis_Script.sh | parallel --jobs "+str(NCores)+" --joblog "+OutFolder+"parallel1_"+OutFilePre+Group+".log\n"
	Line += "date >> "+OutFolder+OutFilePre+Group+".log\n"
	OutScriptGroup.append(Line)
	if RoundNum != NRounds:
		#******tseq_placer_dup.py (merge)
		Line = ScriptFolder+"tseq_placer_dup.py "+LLFileName+" "+ParDictFileName+" "+DataPath+"_exons/ "+OutFilePre+"e"+str(RoundNum)+"_ none same "+OutFilePre+"e"+str(RoundNum)+"_ "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" 0.10 merge >> "+OutFolder+OutFilePre+Group+".log\n"
	else:
		#********tseq_placer_dup.py (separate)
		Line = ScriptFolder+"tseq_placer_dup.py "+LLFileName+" "+OutFolder+OutFilePre+"_round"+str(RoundNum-1)+"/"+OutFilePre+"rnd"+str(RoundNum-1)+"_pardict_new.txt "+DataPath+"_exons/ "+OutFilePre+"e"+str(RoundNum)+"_ none same "+OutFilePre+"e"+str(RoundNum)+"_ "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+OutFolder+OutFilePre+"_round"+str(RoundNum-1)+"/ "+OutFilePre+"rnd"+str(RoundNum-1)+"_ _allseqs_al 0.10 separate >> "+OutFolder+OutFilePre+Group+".log\n"
	Line += "cat "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Analysis_Script.sh | parallel --jobs "+str(NCores)+" --joblog "+OutFolder+"parallel2_"+OutFilePre+Group+".log\n"
	Line += "date >> "+OutFolder+OutFilePre+Group+".log\n"
	OutScriptGroup.append(Line)
	if RoundNum != NRounds:
		#******tcontigs_to_fixed_paralogs.py (intermediate)
		Line = ScriptFolder+"tcontigs_to_fixed_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+DataPath+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ "+DataPath+"_exons/ "+OutFilePre+"e"+str(RoundNum+1)+"_ "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ScriptFolder2+" "+OGFileName+" intermediate none 0 "+OutFilePre+"e"+str(RoundNum)+"_ >> "+OutFolder+OutFilePre+Group+".log\n\n"
	else:
		#********tcontigs_to_fixed_paralogs.py (final/finalpoly)
		if Group in PolyList:
			Line = ScriptFolder+"tcontigs_to_fixed_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+DataPath+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ same none "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ScriptFolder2+" "+OGFileName+" finalpoly "+DataFolder+Group+"/"+GroupDict[Group]+"_poly.txt 2.5 "+OutFilePre+"e"+str(RoundNum)+"_ >> "+OutFolder+OutFilePre+Group+".log\n\n"
		else:
			Line = ScriptFolder+"tcontigs_to_fixed_paralogs.py "+DataPath+"_paralogs/"+OutFilePre+"p"+str(RoundNum)+"_Contig_Groups.txt "+DataPath+"_paralogs/ "+OutFilePre+"p"+str(RoundNum)+"_ "+DataPath+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ same none "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ScriptFolder2+" "+OGFileName+" final none 0 "+OutFilePre+"e"+str(RoundNum)+"_ >> "+OutFolder+OutFilePre+Group+".log\n\n"
	OutScriptGroup.append(Line)
	OutScriptGroupName = OutFolder+OutFilePre+Group+"_"+str(RoundNum)+".sh"
	OutFileWriting(OutScriptGroupName, OutScriptGroup)
	OutLine = "sbatch "+OutScriptGroupName
	SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
	SBatchList.append(SBatchOut.strip('\r').strip('\n').split(" ")[-1])
print("%d separate jobs were submitted to analyze the groups for round %d.\n" % (len(SBatchList), RoundNum))
sys.stderr.write("%d separate jobs were submitted to analyze the groups for round %d.\n" % (len(SBatchList), RoundNum))

#Submitting the next script
OutScript = ['#! /bin/bash\n#SBATCH -J '+OutFilePre+str(RoundNum)+'_comb\n#SBATCH -t 1:00:00\n\n\n']
Line = "rm -r  "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"\nmkdir "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"\n"
OutScript.append(Line)
if RoundNum != NRounds:
	#add the new sequences to the original backbone trees for the first round
	if RoundNum == 1:
		Line = "\n\n\n"+ScriptFolder+"tparalog_combiner.py "+GroupListFileName+" "+OGFileName+" "+DataFolder+" "+DataFolder1Pre+" "+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ "+AlFolder+" "+AlFilePre+" "+AlFilePost+" "+ScriptFolder2+" "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ParDictFileName+" "+str(NCores)+" Array "+OutFolder+OutFilePre+"_"+str(RoundNum+1)+"_group_scripts.sh >> "+OutFolder+OutFilePre+Group+".log\n"
	#and add them to the backbone trees from the previous round for every subsequent round
	else:
		Line = "\n\n\n"+ScriptFolder+"tparalog_combiner.py "+GroupListFileName+" "+OGFileName+" "+DataFolder+" "+DataFolder1Pre+" "+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ "+OutFolder+OutFilePre+"_round"+str(RoundNum-1)+"/ "+OutFilePre+"rnd"+str(RoundNum-1)+"_ _allseqs_al "+ScriptFolder2+" "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ParDictFileName+" "+str(NCores)+" Array "+OutFolder+OutFilePre+"_"+str(RoundNum+1)+"_group_scripts.sh >> "+OutFolder+OutFilePre+Group+".log\n"
	OutScript.append(Line)
	#writing a new script for tparalog_combiner.py to call
	OutScript2 = [ ]
	Line = "#! /bin/bash\n"
	Line += ScriptFolder+"tccm_group_array_subscript.py "+" ".join([ScriptFolder, DataFolder, DataFolder1Pre, DataFolder2Pre, OutFolder, OutFilePre, LLFileName, ParDictFileName, OGFileName, GroupListFileName, PolyListFileName, IndListFileName, AlFolder, AlFilePre, AlFilePost, str(NRounds), str(RoundNum+1)])+"\n"
	OutScript2.append(Line)
	OutFileName = OutFolder+OutFilePre+"_"+str(RoundNum+1)+"_group_scripts.sh"
	OutFileWriting(OutFileName, OutScript2)
else:
	#if it is the last round, we do not want to analyze the trees, so we run tparalog_combiner.py in Parallel mode so that it will produce a script (which we don't run) instead of starting the trees directly.
	Line = "\n\n\n"+ScriptFolder+"tparalog_combiner.py "+GroupListFileName+" "+OGFileName+" "+DataFolder+" "+DataFolder1Pre+" "+DataFolder2Pre+OutFilePre+"_sequences/ "+OutFilePre+"s"+str(RoundNum)+"_ "+OutFolder+OutFilePre+"_round"+str(RoundNum-1)+"/ "+OutFilePre+"rnd"+str(RoundNum-1)+"_ _allseqs_al "+ScriptFolder2+" "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/ "+OutFilePre+"rnd"+str(RoundNum)+"_ "+ParDictFileName+" "+str(NCores)+" Parallel >> "+OutFolder+OutFilePre+Group+".log\n"
	OutScript.append(Line)
	Line = "chmod u+x "+OutFolder+OutFilePre+"_round"+str(RoundNum)+"/"+OutFilePre+"rnd"+str(RoundNum)+"_file_moving_script.sh\n"
	Line += OutFolder+OutFilePre+"_round"+str(RoundNum)+"/"+OutFilePre+"rnd"+str(RoundNum)+"_file_moving_script.sh\n\n\n"
	OutScript.append(Line)
	Line = "sbatch "+OutFolder+OutFilePre+"_contig_classification_script2.sh\n"
	OutScript.append(Line)
OutFileName = OutFolder+OutFilePre+str(RoundNum)+"combined_analysis_script.sh"
OutFileWriting(OutFileName, OutScript)
#this will take place only after the job arrays for all of the individuals have been run.
OutLine = "sbatch -d afterok:"+":".join(SBatchList)+" "+OutFileName
SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
print("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
sys.stderr.write("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
