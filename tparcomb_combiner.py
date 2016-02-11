#! /usr/bin/env python

#tparcomb_combiner.py version 1.0 7 Nov. 2015 Abby Moore
#This script takes the results from multiple runs of tparalog_combiner.py and combines
#them into one final file.  It makes equivalents of all of the output of tparalog_combiner.py:
#Ambigs_per_Seq.txt, Contig_Fates.txt, Ind_Seq_Info.txt, Locus_Paralog_List.txt, 
#perc_ambig.txt, and Seqs_per_Locus.txt.
#When it is combining these files, it somehow needs to deal with having multiple
#sequences per locus.  And it needs to be able to flag the instances of when multiple
#sequences per locus arose in the same round, so they are probably real, and in different rounds,
#which means that they are suspicious and should perhaps be combined.
#It also makes a combined pardict and three types of combined sequence files:
#all sequences for the locus, the best sequences for the locus (all sequences longer
#then 0.75*mean sequence length), and separate files for each (major) paralog.
#Then it writes scripts to add these to the backbone alignments, align them,
#and make trees.
#version 1.1 31 Dec. 2015
#modified so that it doesn't make trees and so that it makes a list of individuals that have
#double sequences for some loci.

#format of Ambigs_per_Seq.txt (has a header line, starting with "Locus" that identifies the groups)
'''
Locus: amk [0]
Paralog: amk1 [1]
Group_#_Ambigs_in_Good_Seqs: 0,2 [2 and subsequent even columns--comma delimitted]
Group_#_Ambigs_in_Bad_Seqs: 14,55,66,17,13,6,22 [3 and subsequent odd columns--comma delimitted]
Columns are blank if there are no sequences for that paralog from that group
'''
#format of Contig_Fates.txt (has a header line, starting with "Locus")
'''
Locus: amk [0]
Group: Anacampserotaceae [1]
Contig: Anacampseros_albidiflora_102-K21-NODE_12_length_535_cov_2.34436_exons [2]
Fate: good_amk1 [3]
'''
#format of Ind_Seq_Info.txt (has a header line, starting with "Locus" that identifies the individuals)
'''
Locus: amk [0]
Paralog: Ambig_amk_none [1]
individual information: 12 seqs, overlap: 630, 59 ambigs, length: 657, using: no [2 and all subsequent columns]
Columns are blank if there are no sequences for that paralog from that group
'''
#format of Locus_Paralog_List.txt (no header)
'''
Locus: amk [0]
Paralog: amk1 [1]
'''
#format of perc_ambig.txt (has a header line, starting with "Locus" that identifies the individuals)
'''
Locus: amk [0]
Paralog: Ambig_amk_none [1]
individual information: n/a or 0.0898 [2 and all subsequent columns]
'''
#format of Seqs_per_Locus.txt (has a header line, starting with "Locus" that identifies the groups)
'''
Locus: amk [0]
Paralog: amk2 [1]
Group_Total_Good_Seqs: 2 [2 and subsequent even columns]
Group_Total_Bad_Seqs: 0	[3 and subsequent odd columns]
'''
#format of Paralog_Info.txt from the tcontig_selection.py (header starting with Locus--right now it starts with Paralog
#and doesn't have a Locus column***)
#Although I guess this may not need to be added to anything, since there should only be one copy.
'''
Locus: amk [0]
Paralog: amk1 [1]
Segment_Number: 2 [2]
Number_of_Contigs: 11 [3]
Segment_Length: 260 [4]
List_of_Groups: Montiaceae [5, comma delimitted]
'''


import sys
from collections import defaultdict
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

Usage = '''
tparcomb_combiner.py combines the information files and sequence files
produced by the various iterations of tparalog_combiner.py and 
tcontig_selection.py and writes scripts to analyze these.  It also flags 
suspicious paralogs, in which different rounds each produced sequences from the
same paralog.
[path to the input files]
[tab-delimitted list with the folder names (tab) prefixes for the input files]
[contig folder]
[prefix for contig files, or "none", if none]
[output folder]
[prefix for output folder]
[folder in which backbone alignments are found]
[prefix for backbone alignments, or "none", if none]
[suffix for backbone alignments, or "none", if none]
[path to folder where scripts are found, or "none", if none]
[pardict file name]
[number of cores for parallelization]
'''

'''
tparcomb_combiner.py InFilePath InFileGroupList IndDictFileName OutGroupDictFN ContigFolder ContigFilePre OutFolder OutFilePre AlFolder AlFilePre AlFilePost ScriptPath ParDictFileName NCores
tparcomb_combiner.py ~/transcriptomes/sandbox/amk/ ~/transcriptomes/sandbox/amk/amk_final/InFileList.txt ~/transcriptomes/general/Ln1_inds.txt ~/transcriptomes/general/outgroup_list_new.txt ~/transcriptomes/sandbox/amk/amk_contigsplit/ amkcs_ ~/transcriptomes/sandbox/amk/amk_final/ amkfi_ ~/transcriptomes/general/phot_optimal_trees/ new_ _al none ~/transcriptomes/general/pardict6.txt 2
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 15:
	sys.exit("ERROR!  This script requires 14 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
InFilePath = sys.argv[1]
if InFilePath[-1] != "/":
	InFilePath += "/"
InFileGroupList = sys.argv[2]
IndDictFileName = sys.argv[3]
OutGroupDictFN = sys.argv[4]
ContigFolder = sys.argv[5]
if ContigFolder[-1] != "/":
	ContigFolder += "/"
ContigFilePre = sys.argv[6]
if ContigFilePre == "none":
	ContigFilePre = ""
OutFolder = sys.argv[7]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[8]
AlFolder = sys.argv[9]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[10]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[11]
if AlFilePost == "none":
	AlFilePost = ""
ScriptPath = sys.argv[12]
if ScriptPath[-1] != "/":
	ScriptPath += "/"
if ScriptPath == "none/":
	ScriptPath = ""
ParDictFileName = sys.argv[13]
NCores = sys.argv[14]

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
	#This is FileGroupDict

#ListDictFromFileList makes a dictionary from multiple tab-delimited files, where the first
#column is the key and the second is the value.  Each key has multiple values, and they
#will be made into a list
#modified from ListDictFromFile in tparalog_combiner.py
def ListDictFromFileList(FileList):
	TempDict = defaultdict(list)
	for InFileName in FileList:
		InFile = open(InFileName, 'rU')
		for Line in InFile:
			Line = Line.strip('\r').strip('\n').split('\t')
			TempDict[Line[0]].append(Line[1])
		InFile.close()
	for Locus in TempDict:
		TempList = TempDict[Locus]
		TempDict[Locus] = list(set(TempList))
	print("%d files were read and combined into a dictionary.\n" % (len(FileList)))
	sys.stderr.write("%d files were read and combined into a dictionary.\n" % (len(FileList)))
	return TempDict
	#This is LPDict

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
	
#DictFromFile makes a dictionary from a tab-delimited file, where the keys are in columns KC, and
#the values are in column VC
#a modified version of the function in from tbaits_intron_removal.py
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
	#This is IndDict

#FileListMaking makes a list of files from a path name, a dictionary of folder name: file prefix,
#whatever goes in between folder name and file prefix (or nothing), whatever goes between the file prefix
#and the file name (or nothing), and the file name.
def FileListMaking(PathName, FilePart1, FilePart2, FileDict, FileName):
	ListTemp = [ ]
	for GroupName in sorted(FileDict.keys()):
		NameTemp = PathName+GroupName+"/"+FilePart1+FileDict[GroupName]+FilePart2+FileName
		ListTemp.append(NameTemp)
	return ListTemp
	#This is various file lists

#This is based on MRScriptWriter from tbaits_introns_removal.py and AMScriptWriter from tclade_finder.py
#It writes a script that adds the sequences to the file of pre-existing
#sequences for that locus and then aligns them.
#from tparalog_combiner.py
def AMRScriptWriter(OLList, OGDict, Folder, Prefix, AFolder, APre, APost, Path):
	OutList1 = ["#! /bin/bash\n"]
	OutList2 = ["#! /bin/bash\n"]
	OutFileName2 = Folder+Prefix+"Locus_Analysis_Subscript.sh"
	#then deal with the alignments
	for Locus in OLList:
		Line = "rm "+Folder+Prefix+Locus+"_allbest_al.fa\n"
		OutList1.append(Line)
		Line = "mafft --localpair --add "+Folder+Prefix+Locus+"_combined_best.fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+Prefix+Locus+"_allbest_al.fa"
		OutList2.append(Line+"\n")
		Line = "rm "+Folder+Prefix+Locus+"_allseqs_al.fa\n"
		OutList1.append(Line)
		Line = "mafft --localpair --add "+Folder+Prefix+Locus+"_combined_all.fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+Prefix+Locus+"_allseqs_al.fa"
		OutList2.append(Line+"\n")
	#writing the script for parallel to analyze
	OutFileWriting(OutFileName2, OutList2)
	Line = "chmod u+x "+OutFileName2+"\n"
	Line += "cat "+OutFileName2+" | parallel --jobs "+NCores+" --progress\n"
	OutList1.append(Line)
	return OutList1
	#This is AnalysisScript.

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

########################################################################################################################################################################


#reading the list of groups
FileGroupDict = NoneDictFromFile(InFileGroupList)
#making a dictionary telling which individual belongs to which group
IndDict = DictFromFile(IndDictFileName, 1, 0)
#making the dictionary telling which sequence is the outgroup for each locus
OutGroupDict = DictFromFile(OutGroupDictFN, 0, 1)
#making the dictionary saying which are the official paralogs for each locus
OrigParalogDict = ListDictFromFile(ParDictFileName, 0, 2)

#Locus_Paralog_List.txt
#making a list of the file names
LPFileList = FileListMaking(InFilePath, "", "", FileGroupDict, "_Locus_Paralog_List.txt")
#combining them into a dictionary
LPDict = ListDictFromFileList(LPFileList)#LPDict[Locus] = list of Paralogs
#making the combined output file
OutFileName = OutFolder+OutFilePre+"Locus_Paralog_List.txt"
OutList = [ ]
for Locus in sorted(LPDict.keys()):
	for Paralog in sorted(LPDict[Locus]):
		Line = Locus+"\t"+Paralog+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#############parsing the Ind_Seq_Info.txt files
#preparing in the dictionary
ISIDict = defaultdict(dict)#ISIDict[Locus][Paralog][Group][Ind][FileGroup]= contig information
SeqsUsingDict = defaultdict(dict)#SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup]='yes'/'no'/'contigs'
IndList = [ ]#The list of all of the individuals that are actually present
GroupList = [ ]#The list of all of the groups that are actually present
GroupIndDict = defaultdict(list)#The dictionary of all of the individuals that are actually present in each group
LastFileGroup = sorted(FileGroupDict.keys())[-1]
for Locus in LPDict:
	for Paralog in LPDict[Locus]:
		ISIDict[Locus][Paralog] = defaultdict(dict)
		SeqsUsingDict[Locus][Paralog] = defaultdict(dict)
for FileGroup in FileGroupDict:
	InFileName = InFilePath+FileGroup+"/"+FileGroupDict[FileGroup]+"_Ind_Seq_Info.txt"
	RoundNum = FileGroup[-1:]
	#first, read the information from the file
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] == "Locus":
			#the line is the list of keys for which individuals are in which positions
			KeyList = Line
			IndRange = range(2,len(KeyList))
			IndList += KeyList[2:]
		else:
			Locus = Line[0]
			Paralog = Line[1]
			for IndPos in IndRange:
				Ind = KeyList[IndPos]
				Group = IndDict[Ind]
				SeqInfo = Line[IndPos]
				if SeqInfo != "":
					SeqInfo = "round "+RoundNum+", "+SeqInfo
					try:
						ISIDict[Locus][Paralog][Group][Ind][FileGroup] = SeqInfo
					except KeyError:
						ISIDict[Locus][Paralog][Group]=defaultdict(dict)
						SeqsUsingDict[Locus][Paralog][Group] = defaultdict(dict)
						ISIDict[Locus][Paralog][Group][Ind][FileGroup] = SeqInfo
						GroupList.append(Group)
					if SeqInfo[-2:] == "no":
						if FileGroup == LastFileGroup:
							SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] = 'contigs'
						else:
							SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] = 'no'
					elif SeqInfo[-3:] == "yes":
						SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] = 'yes'
					else:
						print("ERROR!!!  it is not clear if this sequence has been used or not: %s" % (SeqInfo))
	InFile.close()
#condensing the IndList to come up with the list of individuals that are actually present:
IndList = sorted(list(set(IndList)))
#condensing the GroupList to come up with the list of groups that are actually present:
GroupList = sorted(list(set(GroupList)))
#making the GroupIndDict to show which individuals are actually present for which groups:
for Ind in IndList:
	GroupIndDict[IndDict[Ind]].append(Ind)
for Group in GroupIndDict:
	ListTemp = GroupIndDict[Group]
	ListTemp = sorted(ListTemp)
	GroupIndDict[Group] = ListTemp
#making a list of the FileGroups we're using for each paralog
ParsUsingDict = defaultdict(dict)#ParsUsingDict[Locus][Paralog][Group] = list of FileGroups
RedoList = defaultdict(dict)#RedoList[Locus][Paralog][Group] = list of individuals
for Locus in SeqsUsingDict:
	for Paralog in SeqsUsingDict[Locus]:
		ParsUsingDict[Locus][Paralog]=defaultdict(dict)
		#before I can write the information for the individuals, I need to see how many times that paralog came up
		for Group in GroupList:
			FGListG = [ ]#the list of file groups for which we have sequences we used for that group
			FGListC = [ ]#the list of file groups for which we have contigs we used for that group
			FGListB = [ ]#the list of file groups for which we have either sequences or contigs
			#go through all of the individuals in all of the file groups
			for Ind in SeqsUsingDict[Locus][Paralog][Group]:
				for FileGroup in SeqsUsingDict[Locus][Paralog][Group][Ind]:
					#if we are using the sequence from that individual from that FileGroup,
					if SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] == 'yes':
						#add that FileGroup to the list
						FGListG.append(FileGroup)
						FGListB.append(FileGroup)
					elif SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] == 'contigs':
						FGListC.append(FileGroup)
						FGListB.append(FileGroup)
			FGListG = sorted(list(set(FGListG)))
			ParsUsingDict[Locus][Paralog][Group]['good'] = FGListG
			FGListC = sorted(list(set(FGListC)))
			ParsUsingDict[Locus][Paralog][Group]['contig'] = FGListC
			FGListB = sorted(list(set(FGListB)))
			ParsUsingDict[Locus][Paralog][Group]['both'] = FGListB
#writing the information to the files, while lining up multiple SeqInfos for the same paralog correctly
#ISIDict = defaultdict(dict)#ISIDict[Locus][Paralog][Group][Ind][FileGroup]= contig information
#SeqsUsingDict = defaultdict(dict)#SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup]=True=1/False=0
OutFileName = OutFolder+OutFilePre+"Ind_Seq_Info.txt"
OutList = ["Locus\tParalog\t"+"\t".join(["\t".join(GroupIndDict[Group]) for Group in GroupList])+"\n"]
#also making a dictionary of the overall paralogs we're using for each group
for Locus in sorted(ISIDict.keys()):
	for Paralog in sorted(ISIDict[Locus].keys()):
		Line = Locus+"\t"+Paralog+"\t"
		#before I can write the information for the individuals, I need to see how many times that paralog came up
		for Group in GroupList:
			#if we just have a single FileGroup with good sequences, use it without further complications
			if len(ParsUsingDict[Locus][Paralog][Group]['both']) == 1:
				FGUsing = ParsUsingDict[Locus][Paralog][Group]['both'][0]
				for Ind in GroupIndDict[Group]:
					try:
						Line += ISIDict[Locus][Paralog][Group][Ind][FGUsing]+"\t"
					except KeyError:
						Line += "\t"
			#if there are multiple FileGroups with sequences we're using
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) > 1:
				for Ind in GroupIndDict[Group]:
					#add the information from each FileGroup
					IndGroups = 0
					LinePart = ""
					for FileGroup in ParsUsingDict[Locus][Paralog][Group]['both']:
						try:
							LinePart += "/"+ISIDict[Locus][Paralog][Group][Ind][FileGroup]
							IndGroups += 1
						except KeyError:
							"do nothing"
					#if the individual has more than one group, add them to the RedoList.
					if IndGroups > 1:
						try:
							RedoList[Locus][Paralog][Group].append(Ind)
						except KeyError:
							RedoList[Locus][Paralog] = defaultdict(list)
							RedoList[Locus][Paralog][Group].append(Ind)
					#remove the first slash and add a tab, and add it to the whole line
					Line += LinePart[1:]+"\t"
				#print("Multiple sequences for Group %s for paralog %s.\n" % (Group, Paralog))
			#but, if the group isn't present for that paralog,
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) == 0:
				#insert the correct number of tabs
				Line += "\t"*len(GroupIndDict[Group])
		OutList.append(Line+"\n")
OutFileWriting(OutFileName, OutList)
print("%d loci had paralogs that had multiple sequences in at least one group.  These loci need to be redone.\n" % (len(RedoList.keys())))
sys.stderr.write("%d loci had paralogs that had multiple sequences in at least one group.  These loci need to be redone.\n" % (len(RedoList.keys())))

#printing the RedoList
OutList = ['Locus\tParalog\tGroup\tIndividuals\n']
for Locus in sorted(RedoList.keys()):
	for Paralog in sorted(RedoList[Locus].keys()):
		for Group in sorted(RedoList[Locus][Paralog].keys()):
			Line = Locus+"\t"+Paralog+"\t"+Group+"\t"+(", ").join(RedoList[Locus][Paralog][Group])+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"RedoList.txt"
OutFileWriting(OutFileName, OutList)

##I need to figure out which sequence names belong to divided sequences
PrunedLPDict = defaultdict(dict)
#dictionary of ambiguous paralogs, so I can figure out which groups have them.
AmbigParDict = defaultdict(dict)
for Locus in ISIDict:
	PrunedLPDict[Locus] = defaultdict(list)
	for Paralog in ISIDict[Locus]:
		#There are three possibilities for the paralogs:
		#First, they can be ambiguous, in which case they won't be classified in a particular paralog.
		if Paralog[:5] == "Ambig":
			PrunedLPDict[Locus][Paralog] = [Paralog]
			AmbigParDict[Locus][Paralog] = [ ]
		#Second, they can be the original paralog, in which case they are classified to the original paralog.
		elif (Paralog in OrigParalogDict[Locus]) == True:
			PrunedLPDict[Locus][Paralog].append(Paralog)
		#Third, the original paralog can have been divided in that group.
		else:
			ParalogTemp = Paralog.split("_")
			#first, trying the first part of the name
			if (len(ParalogTemp) == 2) and (ParalogTemp[0] in OrigParalogDict[Locus]):
				PrunedLPDict[Locus][ParalogTemp[0]].append(Paralog)
			#then seeing if the original paralog name could have two parts, separated by an underscore
			elif (len(ParalogTemp) == 3) and (ParalogTemp[0]+"_"+ParalogTemp[1] in OrigParalogDict[Locus]):
				PrunedLPDict[Locus][ParalogTemp[0]+"_"+ParalogTemp[1]].append(Paralog)
			#then, assuming that the paralog underwent multiple rounds of separation (so that this is part 1 of part 1, for example)
			elif (len(ParalogTemp) > 2) and (ParalogTemp[0] in OrigParalogDict[Locus]):
				PrunedLPDict[Locus][ParalogTemp[0]].append(Paralog)
			elif (len(ParalogTemp) > 2) and (ParalogTemp[0]+"_"+ParalogTemp[1] in OrigParalogDict[Locus]):
				PrunedLPDict[Locus][ParalogTemp[0]+"_"+ParalogTemp[1]].append(Paralog)
			else:
				print("ERROR!!  This paralog does not fit the pattern!! %s." % (Paralog))
				sys.stderr.write("ERROR!!  This paralog does not fit the pattern!! %s." % (Paralog))
#and writing a Seq_Fates.txt file with the summary of the various lists
OutFileName1 = OutFolder+OutFilePre+"Seq_Fates_All.txt"
OutFileName2 = OutFolder+OutFilePre+"Seq_Fates_wholeseqs.txt"
OutFileName3 = OutFolder+OutFilePre+"Seq_Fates_contigs.txt"
OutList1 = ["Locus\tParalog\t"+"\t".join(["\t".join(GroupIndDict[Group]) for Group in GroupList])+"\n"]
OutList2 = ["Locus\tParalog\t"+"\t".join(["\t".join(GroupIndDict[Group]) for Group in GroupList])+"\n"]
OutList3 = ["Locus\tParalog\t"+"\t".join(["\t".join(GroupIndDict[Group]) for Group in GroupList])+"\n"]
#also making a dictionary of the overall paralogs we're using for each group
for Locus in sorted(PrunedLPDict.keys()):
	for Paralog in sorted(PrunedLPDict[Locus].keys()):
		Line1 = Locus+"\t"+Paralog+"\t"
		Line2 = Locus+"\t"+Paralog+"\t"
		Line3 = Locus+"\t"+Paralog+"\t"
		for Group in GroupList:
			for Ind in GroupIndDict[Group]:
				NumYes = 0
				NumContigs = 0
				for SubParalogName in PrunedLPDict[Locus][Paralog]:
					try:
						for FileGroup in SeqsUsingDict[Locus][SubParalogName][Group][Ind]:
							if SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] == 'yes':
								NumYes += 1
							elif SeqsUsingDict[Locus][Paralog][Group][Ind][FileGroup] == 'contigs':
								NumContigs += 1
					except KeyError:
						"do nothing"
				Line1 += str(NumYes)+"_"+str(NumContigs)+"\t"
				Line2 += str(NumYes)+"\t"
				Line3 += str(NumContigs)+"\t"
		OutList1.append(Line1[:-1]+"\n")
		OutList2.append(Line2[:-1]+"\n")
		OutList3.append(Line3[:-1]+"\n")
OutFileWriting(OutFileName1, OutList1)
OutFileWriting(OutFileName2, OutList2)
OutFileWriting(OutFileName3, OutList3)

#filling out and writing the AmbigParDict
OutList = [ ]
for Locus in AmbigParDict:
	for Paralog in AmbigParDict[Locus]:
		Line = Locus+"\t"+Paralog+"\t"+",".join(ISIDict[Locus][Paralog].keys())+"\n"
		OutList.append(Line)
		Line = "\t".join([Locus, Paralog, "All"])+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Ambig_Paralog_List.txt"
OutFileWriting(OutFileName, OutList)


##############Ambigs_per_Seq.txt
ApSDict = defaultdict(dict)#ApSDict[Locus][Paralog][Group][FileNum]['Good'/'Bad'] = comma-delimitted list
for Locus in LPDict.keys():
	ApSDict[Locus] = defaultdict(dict)
for FileGroup in FileGroupDict:
	InFileName = InFilePath+FileGroup+"/"+FileGroupDict[FileGroup]+"_Ambigs_per_Seq.txt"
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] == "Locus":
			KeyList = [Item.split("_#")[0] for Item in Line]
			LinePosList = range(2,len(KeyList), 2)
		else:
			Locus = Line[0]
			Paralog = Line[1]
			for LinePos in LinePosList:
				Group = KeyList[LinePos]
				Good = Line[LinePos]
				Bad = Line[LinePos+1]
				if (Good != "") or (Bad != ""):
					try:
						ApSDict[Locus][Paralog][Group][FileGroup]['Good'] = Good
						ApSDict[Locus][Paralog][Group][FileGroup]['Bad'] = Bad
					except KeyError:
						ApSDict[Locus][Paralog][Group] = defaultdict(dict)
						ApSDict[Locus][Paralog][Group][FileGroup]['Good'] = Good
						ApSDict[Locus][Paralog][Group][FileGroup]['Bad'] = Bad
	InFile.close()
#printing the file
OutFileName = OutFolder + OutFilePre + "Ambigs_per_Seq.txt"
OutList = ["Locus\tParalog\t"+"\t".join([Group+"_#_Ambigs_in_Good_Seqs\t"+Group+"_#_Ambigs_in_Bad_Seqs" for Group in GroupList])+"\n"]
for Locus in sorted(ApSDict.keys()):
	for Paralog in sorted(ApSDict[Locus].keys()):
		Line = Locus+"\t"+Paralog+"\t"
		for Group in GroupList:
			#print("%s: %s: %s\n" % (Locus, Paralog, Group))
			#if we just have a single FileGroup with good sequences, use it without further complications
			if len(ParsUsingDict[Locus][Paralog][Group]['both']) == 1:
				FGUsing = ParsUsingDict[Locus][Paralog][Group]['both'][0]
				Line += ApSDict[Locus][Paralog][Group][FGUsing]['Good']+"\t"+ApSDict[Locus][Paralog][Group][FGUsing]['Bad']+"\t"
			#if there are multiple FileGroups with sequences we're using
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) > 1:
				#add the information from each FileGroup
				GoodList = [ ]
				BadList = [ ]
				for FileGroup in ParsUsingDict[Locus][Paralog][Group]['both']:
					if ApSDict[Locus][Paralog][Group][FileGroup]['Good'] != "":
						GoodList.append(ApSDict[Locus][Paralog][Group][FileGroup]['Good'])
					else:
						GoodList.append("-")
					if ApSDict[Locus][Paralog][Group][FileGroup]['Bad'] != "":
						BadList.append(ApSDict[Locus][Paralog][Group][FileGroup]['Bad'])
					else:
						BadList.append("-")
				Line += "/".join(GoodList)+"\t"+"/".join(BadList)+"\t"
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) == 0:
				#insert the correct number of tabs
				Line += "\t"+"\t"
		Line = Line[:-1]+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#################parsing the Contig_Fates.txt files
#setting up the dictionary I will fill out
CFDict = defaultdict(dict)#CFDict[Locus][Group][FileGroup][Contig] = Fate
for Locus in LPDict.keys():
	for Group in GroupList:
		CFDict[Locus][Group] = defaultdict(dict)
#filling out the dictionary by reading the files:
for FileGroup in FileGroupDict:
	InFileName = InFilePath+FileGroup+"/"+FileGroupDict[FileGroup]+"_Contig_Fates.txt"
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] != "Locus":
			Locus = Line[0]
			Group = Line[1]
			Contig = Line[2]
			Fate = Line[3]
			CFDict[Locus][Group][Contig][FileGroup] = Fate
	InFile.close()
#writing this information to the combined file:
OutFileName = OutFolder+OutFilePre+"Contig_Fates.txt"
OutList = ['Locus\tGroup\tContig\tContig_Fate\n']
for Locus in sorted(CFDict.keys()):
	for Group in sorted(CFDict[Locus].keys()):
		for Contig in sorted(CFDict[Locus][Group].keys()):
			FGUsing = sorted(CFDict[Locus][Group][Contig].keys())[-1]
			Line = Locus+"\t"+Group+"\t"+Contig+"\t"+CFDict[Locus][Group][Contig][FGUsing]+"\n"
			OutList.append(Line)
OutFileWriting(OutFileName, OutList)


#######################dealing with perc_ambig.txt
#preparing in the dictionary
PADict = defaultdict(dict)#PADict[Locus][Paralog][Ind][FileGroup]= percent ambiguities or n/a
for Locus in LPDict.keys():
	for Paralog in LPDict[Locus]:
		PADict[Locus][Paralog] = defaultdict(dict)
#reading the files and writing them to the dictionary
IndList = [ ]
for FileGroup in FileGroupDict:
	InFileName = InFilePath+FileGroup+"/"+FileGroupDict[FileGroup]+"_perc_ambig.txt"
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		#use the header line to read the individuals
		if Line[0] == "Locus":
			KeyList = Line
			IndList += KeyList[2:]
			PARange = range(2,len(KeyList))
		else:
			Locus = Line[0]
			Paralog = Line[1]
			for PANum in PARange:
				Ind = KeyList[PANum]
				PADict[Locus][Paralog][Ind][FileGroup] = Line[PANum]
	InFile.close()
#writing the output
OutFileName = OutFolder+OutFilePre+"perc_ambig.txt"
IndList = sorted(list(set(IndList)))
OutList = ["Locus\tParalog\t"+"\t".join(["\t".join(GroupIndDict[Group]) for Group in GroupList])+"\n"]
for Locus in PADict:
	for Paralog in PADict[Locus]:
		Line = Locus+"\t"+Paralog+"\t"
		for Group in GroupList:
			#if we just have a single FileGroup with good sequences, use it without further complications
			if len(ParsUsingDict[Locus][Paralog][Group]['both']) == 1:
				FGUsing = ParsUsingDict[Locus][Paralog][Group]['both'][0]
				for Ind in GroupIndDict[Group]:
					try:
						Line += PADict[Locus][Paralog][Ind][FGUsing]+"\t"
					except KeyError:
						Line += "n/a\t"
			#if there are multiple FileGroups with sequences we're using, taking the biggest number, because it is likely
			#that the smaller number(s) is(are) due to single sequences that were only properly classified in later rounds,
			#so the larger number of ambiguities is the number for the main sequence.
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) > 1:
				for Ind in GroupIndDict[Group]:
					IndAmbigList = [ ]
					#add the information from each FileGroup
					for FileGroup in ParsUsingDict[Locus][Paralog][Group]['both']:
						try:
							IndAmbigList.append(PADict[Locus][Paralog][Ind][FileGroup])
						except KeyError:
							"do nothing"
					if len(IndAmbigList) == 0:
						Line += "n/t\t"
					elif len(IndAmbigList) == 1:
						Line += IndAmbigList[0]+"\t"
					elif len(IndAmbigList) > 1:
						IndAmbigList = sorted(IndAmbigList)
						Line += IndAmbigList[-1]+"\t"
			#but, if the group isn't present for that paralog,
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) == 0:
				#insert the correct number of tabs
				Line += "n/a\t"*len(GroupIndDict[Group])
		Line = Line[:-1]+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

################parsing Seqs_per_Locus

#dealing with Seqs_per_Locus.txt
#preparing the dictionary
SpLDict = defaultdict(dict)#SpLDict[Locus][Paralog][Group][FileGroup]['Good'/'Bad'] # good/bad sequences
for Locus in LPDict.keys():
	for Paralog in LPDict[Locus]:
		SpLDict[Locus][Paralog]= defaultdict(dict)
#reading the files and writing them to the dictionary
IndList = [ ]
for FileGroup in FileGroupDict:
	InFileName = InFilePath+FileGroup+"/"+FileGroupDict[FileGroup]+"_Seqs_per_Locus.txt"
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] == "Locus":
			KeyList = [Item.split("_Total_")[0] for Item in Line]
			LinePosList = range(2,len(KeyList), 2)
		else:
			Locus = Line[0]
			Paralog = Line[1]
			for LinePos in LinePosList:
				Group = KeyList[LinePos]
				Good = Line[LinePos]
				Bad = Line[LinePos+1]
				if (Good != "") or (Bad != ""):
					try:
						SpLDict[Locus][Paralog][Group][FileGroup]['Good'] = Good
						SpLDict[Locus][Paralog][Group][FileGroup]['Bad'] = Bad
					except KeyError:
						SpLDict[Locus][Paralog][Group][FileGroup] = defaultdict(dict)
						SpLDict[Locus][Paralog][Group][FileGroup]['Good'] = Good
						SpLDict[Locus][Paralog][Group][FileGroup]['Bad'] = Bad
	InFile.close()
#writing the output
OutFileName = OutFolder + OutFilePre + "Seqs_per_Locus.txt"
OutList = ["Locus\tParalog\t"+"\t".join([Group+"_Total_Good_Seqs\t"+Group+"_Total_Bad_Seqs" for Group in GroupList])+"\n"]
for Locus in sorted(SpLDict.keys()):
	for Paralog in sorted(SpLDict[Locus].keys()):
		Line = Locus+"\t"+Paralog+"\t"
		for Group in GroupList:
			#if we just have a single FileGroup with good sequences, use it without further complications
			if len(ParsUsingDict[Locus][Paralog][Group]['both']) == 1:
				FGUsing = ParsUsingDict[Locus][Paralog][Group]['both'][0]
				Line += SpLDict[Locus][Paralog][Group][FGUsing]['Good']+"\t"+SpLDict[Locus][Paralog][Group][FGUsing]['Bad']+"\t"
			#if there are multiple FileGroups with sequences we're using
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) > 1:
				#sum the good sequences and bad sequences, because, if there are multiple good/bad sequences, they will
				#all have been included
				GoodList = [ ]
				BadList = [ ]
				for FileGroup in ParsUsingDict[Locus][Paralog][Group]['both']:
					if SpLDict[Locus][Paralog][Group][FileGroup]['Good'] != "":
						GoodList.append(SpLDict[Locus][Paralog][Group][FileGroup]['Good'])
					if SpLDict[Locus][Paralog][Group][FileGroup]['Bad'] != "":
						BadList.append(SpLDict[Locus][Paralog][Group][FileGroup]['Bad'])
				if len(GoodList) == 0:
					GoodSeqs = str(0)
				elif len(GoodList) == 1:
					GoodSeqs = GoodList[0]
				elif len(GoodList) > 1:
					GoodSeqs = str(sum([int(Num) for Num in GoodList]))
				if len(BadList) == 0:
					BadSeqs = str(0)
				elif len(BadList) == 1:
					BadSeqs = BadList[0]
				elif len(BadList) > 1:
					BadSeqs = str(sum([int(Num) for Num in BadList]))
				Line += GoodSeqs+"\t"+BadSeqs+"\t"
			elif len(ParsUsingDict[Locus][Paralog][Group]['both']) == 0:
				#insert the correct number of tabs
				Line += "0\t0\t"
		Line = Line[:-1]+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#combining the sequences, first the combined sequences for each locus
NumBestSeqFiles = 0
NumAllSeqFiles = 0
NumParalogFiles = 0
for Locus in LPDict.keys():
	#first the best sequences, because the contigs won't be added here
	OutFileName2 = OutFolder+OutFilePre+Locus+"_combined_best.fa"
	NumBestSeqFiles += 1
	OutFile2 = open(OutFileName2, 'w')
	InFileList2 = FileListMaking(InFilePath, "", "", FileGroupDict, "_"+Locus+"_combined_best.fa")
	for InFileName in InFileList2:
		try:
			InFile = open(InFileName, 'rU')
			for record in SeqIO.parse(InFile, "fasta"):
				SeqIO.write(record, OutFile2, "fasta")
			InFile.close()
		except IOError:
			"no sequence files for this round"
	OutFile2.close()
	#all sequences
	OutFileName1 = OutFolder+OutFilePre+Locus+"_combined_all.fa"
	NumAllSeqFiles += 1
	OutFile1 = open(OutFileName1, 'w')
	#add the full sequences and leave the file open to add the contigs
	InFileList1 = FileListMaking(InFilePath, "", "", FileGroupDict, "_"+Locus+"_combined_all.fa")
	for InFileName in InFileList1:
		try:
			InFile = open(InFileName, 'rU')
			for record in SeqIO.parse(InFile, "fasta"):
				SeqIO.write(record, OutFile1, "fasta")
			InFile.close()
		except IOError:
			"no sequence files for this round"
	#then deal with the separate paralogs
	for Paralog in LPDict[Locus]:
		OutFileName3 = OutFolder+OutFilePre+Locus+"_"+Paralog+".fa"
		NumParalogFiles += 1
		OutFile3 = open(OutFileName3, 'w')
		#first add the good paralog sequences to the file for that paralog
		InFileList3 = FileListMaking(InFilePath, "", "", FileGroupDict, "_"+Locus+"_"+Paralog+".fa")
		for InFileName in InFileList3:
			try:
				InFile = open(InFileName, 'rU')
				for record in SeqIO.parse(InFile, "fasta"):
					SeqIO.write(record, OutFile3, "fasta")
				InFile.close()
			except IOError:
				"no sequence files for this round"
		#then look through the separated contigs
		InFileName = ContigFolder+ContigFilePre+Paralog+".fa"
		try:
			InFile = open(InFileName, 'rU')
			for record in SeqIO.parse(InFile, "fasta"):
				#add them to the file for that paralog
				SeqIO.write(record, OutFile3, "fasta")
				#and to the file with all sequences for that locus
				SeqIO.write(record, OutFile1, "fasta")
			InFile.close()
		except IOError:
			"no sequence files for this record"
		OutFile3.close()
	OutFile1.close()
#output stuff about the files
print("%d files containing the best sequences from each locus were written, with names such as %s (for potential use as backbone alignments).\n" % (NumBestSeqFiles, OutFileName2))
sys.stderr.write("%d files containing the best sequences from each locus were written, with names such as %s (for potential use as backbone alignments).\n" % (NumBestSeqFiles, OutFileName2))
print("%d files containing all of the sequences from each locus were written, with names such as %s (for the final phylogenetic analyses).\n" % (NumAllSeqFiles, OutFileName1))
sys.stderr.write("%d files containing all of the sequences from each locus were written, with names such as %s (for the final phylogenetic analyses).\n" % (NumAllSeqFiles, OutFileName1))
print("%d separate files for each paralog were written, with names such as %s.\n" % (NumParalogFiles, OutFileName3))
sys.stderr.write("%d separate files for each paralog were written, with names such as %s.\n" % (NumParalogFiles, OutFileName3))

####Now to write a script to analyze these files.
AnalysisScript = AMRScriptWriter(LPDict.keys(), OutGroupDict, OutFolder, OutFilePre, AlFolder, AlFilePre, AlFilePost, ScriptPath)
#concatenating the pardicts
PDFileList = FileListMaking(InFilePath, "", "", FileGroupDict, "_pardict.txt")
Line = "cat "+" ".join(PDFileList)+" "+ParDictFileName+" > "+OutFolder+OutFilePre+"pardict_new.txt\n"
AnalysisScript.append(Line)
OutFileName = OutFolder+OutFilePre+"Locus_Analysis_Script.sh"
OutFileWriting(OutFileName, AnalysisScript)
#The paralogs do not need to be aligned, because that will happen with tparcomb_final.py
