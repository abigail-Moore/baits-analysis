#! /usr/bin/env python

#tfig_making_prep.py version 1.0 6 Feb. 2017
#This script takes various files from other scripts and combines them to make figures.

#input files:
#c2p1gt2_seq_length.txt, produced by tnotung_homolog_parsing.py
#loci are columns, individuals are rows, values are sequence length per individual per locus
#Of course, this is before any pruning has taken place.

#c2p1gt2_dup_pos_dict.txt, produced by tnotung_homolog_parsing.py
#first column is gene family name, second column is duplication name, third column is the individuals in the duplication, fourth column is the sister group

#c2p1gt2_Inds_per_Paralog.txt, produced by tparalog_selector.py
#first column is group name, remaining columns are the number of individuals with sequences for each locus

import sys
from collections import defaultdict
import numpy

'''
tfig_making_prep.py InFolder InFilePre OutFolder IndListFN LocusListListFN IndGroupListFN
'''

Usage = '''
tfig_making_prep.py makes figures from the output of tnotung_homolog_parsing.py
and tparalog_selector.py
[folder in which the output from these files is found]
[prefix for input and output files files]
[output folder]
[list of individuals to include]
[list of files with locus lists]
[list of individuals and the groups they are in--can include more individuals 
than are in the final list of individuals]
'''

if len(sys.argv) < 7:
	sys.exit("ERROR!!!  This script requires 6 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
if InFolder[-1] != "/":
	InFolder += "/"
InFilePre = sys.argv[2]
if InFilePre == "none":
	InFilePre = ""
OutFolder = sys.argv[3]
if OutFolder == "same":
	OutFolder = InFolder
elif OutFolder[-1] != "/":
	OutFolder += "/"
IndListFN = sys.argv[4]
LocusListListFN = sys.argv[5]
IndGroupListFN = sys.argv[6]

print(" ".join(sys.argv))

#################################################################################

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
#from tseq_placer_dup.py
def CaptureColumn(FileName, ColNum):
	TempList = [ ]
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r').split('\t')
		TempList.append(Line[ColNum])
	InFile.close()
	#print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	#sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is LocusList

#DictFromFile makes a dictionary from a tab-delimited file, where the keys are in columns KC, and
#the values are in column VC
#a modified version of the function in from tbaits_intron_removal.py
#from tcontig_selection.py
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

#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
#from tbaits_intron_removal.py
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	print("Output file %s written.\n" % (FileName))
	#sys.stderr.write("Output file %s written.\n" % (FileName))

#################################################################################

#the various dictionaries
SeqLenDict = defaultdict(dict)#SeqLenDict[IndName][LocName] = LocLength--from InFilePreseq_length.txt
DupNumDict = defaultdict(int)#DupNumDict[GenFam] = NumDups--from InFilePredup_pos_dict.txt
IndpParDict = defaultdict(dict)#IndpParDict[Locus][GroupName] = NumInds--from InFilePreInds_per_Paralog.txt

#getting the list of individuals we want to look at
IndList = CaptureColumn(IndListFN, 0)

#getting the list of groups in which those individuals are found
IndGroupDict = DictFromFile(IndGroupListFN, 1, 0)

#reading the file showing the locus length for each individual
InFileName = InFolder+InFilePre+"seq_length.txt"
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	if Line[0] == "IndividualName":
		KeyList = Line
		KeyRange = range(1,len(KeyList))
	else: 
		IndName = Line[0]
		for KeyNum in KeyRange:
			SeqLenDict[IndName][KeyList[KeyNum]] = Line[KeyNum]
InFile.close()

LocusList = sorted(KeyList[1:])

#calculate mean locus length and number of individuals
LocusInfoDict = defaultdict(dict)
MissingIndList = []
#making a list of all groups (not the groups in one locus)
GroupList = [ ]
for Locus in LocusList:
	NumInds = 0
	LocusLengthList = [ ]
	GroupListTemp = [ ]
	for IndName in IndList:
		try: LLength = SeqLenDict[IndName][Locus]
		#if there is an individual in the locus list that is not in the sequences
		except KeyError:
			MissingIndList.append(IndName)
		if LLength != "":
			NumInds += 1
			LocusLengthList.append(int(LLength))
			GroupListTemp.append(IndGroupDict[IndName])
	LocusInfoDict[Locus]['NumInds'] = NumInds
	LocusInfoDict[Locus]['MeanLen'] = numpy.mean(LocusLengthList)
	GroupListTemp = list(set(GroupListTemp))
	GroupList += GroupListTemp
	LocusInfoDict[Locus]['NumGroups'] = len(GroupListTemp)
	LocusInfoDict[Locus]['LocusGroupList'] = GroupListTemp
#dealing with any missing individuals, so they don't mess things up later
if MissingIndList != [ ]:
	MissingIndList = sorted(list(set(MissingIndList)))
	for IndName in MissingIndList:
		print ("%s is not found in the locus file.\n" % (IndName))
		IndList.remove(IndName)
		del SeqLenDict[IndName]
#making a better list of all groups so that it can be used in the output
GroupList = sorted(list(set(GroupList)))
	

#printing this information to a file:
#c2p1gt2_length_vs_numinds_all.csv: IndividualName [should be locus_name],mean_length,num_inds
OutFileName = OutFolder+InFilePre+"length_vs_numinds_all.csv"
OutList = ["locus_name,mean_length,num_inds\n"]
for Locus in LocusList:
	Line = ("%s,%.2f,%d\n" % (Locus, LocusInfoDict[Locus]['MeanLen'], LocusInfoDict[Locus]['NumInds']))
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#######reading c2p1gt2_dup_pos_dict.txt and producing c2p1gt2_dups.csv
#c2p1gt2_dup_pos_dict.txt: Locus, Duplication_Name, Duplicated_Individuals, Sister_Group
InFileName = InFolder+InFilePre+"dup_pos_dict.txt"
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	if Line[0] != "Locus":
		DupNumDict[Line[0]] += 1
InFile.close()

GenFamList = sorted(DupNumDict.keys())

#printing this information to a file:
#c2p1gt2_dups.csv: Locus [gene family], Num_Dups, Num_Copies
OutFileName = OutFolder+InFilePre+"dups.csv"
OutList = ['Gene_Fam,Num_Dups,Num_Copies\n']
for GenFam in GenFamList:
	Line = ("%s,%d,%d\n" % (GenFam, DupNumDict[GenFam], DupNumDict[GenFam]+1))
	OutList += Line
OutFileWriting(OutFileName, OutList)

########c2p1gt2_Inds_per_Paralog.txt and producing c2p1gt2_Inds_per_Paralog.csv
#c2p1gt2_Inds_per_Paralog.txt: Group, columns for the various loci
InFileName = InFolder+InFilePre+"Inds_per_Paralog.txt"
GroupKeyList = [ ]
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	if Line[0] == "Group":
		KeyList = Line
		KeyRange = range(1,len(KeyList))
	else: 
		GroupName = Line[0]
		GroupKeyList.append(GroupName)
		for KeyNum in KeyRange:
			IndpParDict[KeyList[KeyNum]][GroupName] = Line[KeyNum]
InFile.close()

#printing this information to a file:
#c2p1gt2_Inds_per_Paralog.csv: Group [meaning locus_name], Anacampserotaceae_st, Basellaceae_st, ..., Total_Sequences_in_All_Groups, Number_of_Groups_with_Sequences
OutFileName = OutFolder+InFilePre+"Inds_per_Paralog.csv"
OutList = ['Group,'+','.join(GroupKeyList)+'\n']
for Locus in LocusList:
	Line = Locus+','+','.join([IndpParDict[Locus][GroupName] for GroupName in GroupKeyList])+'\n'
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)

######making a file with the mean number of loci per individual and the mean number of missing loci
LocusListList = CaptureColumn(LocusListListFN,0)
LocusListList = sorted(LocusListList)
SubsetInfoDict = defaultdict(dict)#SubsetInfoDict[LocusListFN]['NumLoci'/'MeanNumLoci'/'PercMissing']
LocusListDict = defaultdict(dict)#LocusListDict[LocusListFN] = [list of loci in the locus list]
for LocusListFN in LocusListList:
	LocusListDict[LocusListFN] = CaptureColumn(InFolder+LocusListFN, 0)
	NumLoci = len(LocusListDict[LocusListFN])
	SubsetInfoDict[LocusListFN]['NumLoci'] = NumLoci
	TotalSeqs = 0
	for IndName in SeqLenDict:
		#And also getting information about how many loci each individual has
		SubsetInfoDict[LocusListFN][IndName] = 0
		for LocName in LocusListDict[LocusListFN]:
			if SeqLenDict[IndName][LocName] != '':
				TotalSeqs += 1
				SubsetInfoDict[LocusListFN][IndName] += 1
	MeanNumLoci = float(TotalSeqs)/len(IndList)
	SubsetInfoDict[LocusListFN]['MeanNumLoci'] = MeanNumLoci
	PercMissing = 1.0-MeanNumLoci/NumLoci
	SubsetInfoDict[LocusListFN]['PercMissing'] = PercMissing

#printing the summary information to a file
OutFileName = OutFolder+InFilePre+"Locus_Subset_Info.txt"
OutList = ["Subset_Name\tNumber_of_Loci\tMean_Number_of_Loci_per_Individual\tProportion_Missing\n"]
for LocusListFN in LocusListList:
	Line = ("%s\t%d\t%.2f\t%.4f\n" % (LocusListFN, SubsetInfoDict[LocusListFN]['NumLoci'], SubsetInfoDict[LocusListFN]['MeanNumLoci'], SubsetInfoDict[LocusListFN]['PercMissing']))
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#and making a separate file for the information for each individual
OutFileName = OutFolder+InFilePre+"Locus_Subset_Info_Inds.txt"
OutList = ["IndName\t"+"\t".join(LocusListList)+"\n"]
NewIndList = sorted(SeqLenDict.keys())
for IndName in NewIndList:
	Line = IndName
	for LocusListFN in LocusListList:
		Line += "\t"+str(SubsetInfoDict[LocusListFN][IndName])
	OutList.append(Line+"\n")
OutFileWriting(OutFileName, OutList)

##making a file with a bunch of this information combined
#locus, ave_length, num_individs, num_groups, g2, g5, g9, i36, i57, and the various groups (just 0/1 or yes/no for these last ones)
OutFileName = OutFolder+InFilePre+"coverage_data.csv"
OutList = ["Locus,Mean_Length,Num_Inds,Num_Groups,"+",".join(LocusListList)+","+",".join(GroupList)+"\n"]
for Locus in LocusList:
	if LocusInfoDict[Locus]['NumInds'] > 0:
		Line = ("%s,%.2f,%d,%d" % (Locus, LocusInfoDict[Locus]['MeanLen'], LocusInfoDict[Locus]['NumInds'], LocusInfoDict[Locus]['NumGroups']))
		for ListName in LocusListList:
			if Locus in LocusListDict[ListName]:
				Line += ",yes"
			else:
				Line += ",no"
		for GroupName in GroupList:
			if GroupName in LocusInfoDict[Locus]['LocusGroupList']:
				Line += ",yes"
			else:
				Line += ",no"
		OutList.append(Line+"\n")
OutFileWriting(OutFileName, OutList)