#! /usr/bin/env python

#tparalog_combiner.py version 1.1 23 Feb. 2016 Abby Moore
#This script combines the sequences from the various paralogs found by tcontigs_to_paralogs.py
#into a single alignment, along with the backbone sequences.  It also labels the
#sequences according to their paralog.
#version 1.0 23 June 2015
#version 1.1 23 Feb. 2016 modified so that the sequences are ordered in parallel
#**************This changed version has not been tested.***********


#It also combines the following files with information about the paralogs:
#OutFolder+OutFilePre_Ind_Seq_Info.txt (tab-delimitted)
'''
[0] Locus: alaAT
[1] Paralog:Ambig_alaAT2_alaAT1
[2] IndName: Alluaudia_procera_53
[3] #_Sequences: 4
[4] #_Ambiguities:3 
[5] #Seq_Length: 507
[6] Using_Seq: no
'''
#OutFolder+OutFilePre_Paralog_Info.txt (tab-delimitted)
'''
[0] Locus: ppc
[1] Paralog: ppc-1E2
[2] #_Ambiguities_in_Good_Seqs: 0,0,0,0,0,0,0,0,0 (comma-delimitted, can be blank if there are no good sequences)
[3] #_Ambiguities_in_Bad_Seqs: 13,21 (comma-delimitted, can be blank if there are no bad sequences)
[4] #_Good_Seqs: 9
[5] #_Bad_Seqs: 2
'''

from collections import defaultdict
import sys
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
import numpy #to calculate means
import subprocess #We want to be able to talk to the command line.

Usage = '''
tparalog_combiner.py labels sequences according to paralog and combines them.
tparalog_combiner.py
[text file with the names of the different folders, with the name of the sequences
folder in the first column and the shortened name for the subfolders in the 2nd]
[tab-delimitted file with locus in the first column and outgroup sequence name 
in the second column]
[directory in which all of these files are found]
[anything between the name of the main folder and the name of the subfolder or 
"none", if nothing]
[anything after the name of the subfolder, or "none", if nothing--this can just
be "/", if there is no suffix, but there is still a folder]
[prefix for the sequence files, or "none" if none]
[folder of the backbone alignment files]
[prefix of the backbone alignment files, or "none" if none]
[suffix of the backbone alignment files, or "none" if none]
[path to the folder containing scripts, or "none" if none]
[output folder]
[prefix for output files, or "none" if none]
[the name of the pardict file, showing which sequences belong to which paralogs]
[number of cores for script]
[mode: Parallel or Array]
[file name for the next script, if mode is Array]
'''

#tparalog_combiner.py GFileName OGFileName SeqFilePath SeqFolderPre SeqFolderPost SeqFilePre AlFolder AlFilePre AlFilePost ScriptPath OutFolder OutFilePre PFileName NCores Mode[Parallel, Array] NextScript
#tparalog_combiner.py ~/transcriptomes/general/Group_List_sand1.txt ~/transcriptomes/general/outgroup_list_new.txt ~/transcriptomes/sandbox/ none none s1_ ~/transcriptomes/general/phot_optimal_trees/ new_ _al none ~/transcriptomes/sandbox/amk/combined/ amktr1_ ~/transcriptomes/general/pardict6.txt

print("%s\n" % (" ".join(sys.argv)))

ModeList = ['Parallel', 'Array']

if len(sys.argv) < 16:
	sys.exit("ERROR!  This script requires 15 or 16 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
GFileName = sys.argv[1]
OGFileName = sys.argv[2]
SeqFilePath = sys.argv[3]
if SeqFilePath[-1] != "/":
	SeqFilePath += "/"
SeqFolderPre = sys.argv[4]
if SeqFolderPre == "none":
	SeqFolderPre = ""
SeqFolderPost = sys.argv[5]
if SeqFolderPost[-1] != "/":
	SeqFolderPost += "/"
if SeqFolderPost == "none/":
	SeqFolderPost = ""
SeqFilePre = sys.argv[6]
if SeqFilePre == "none":
	SeqFilePre = ""
AlFolder = sys.argv[7]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[8]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[9]
if AlFilePost == "none":
	AlFilePost = ""
ScriptPath = sys.argv[10]
if ScriptPath[-1] != "/":
	ScriptPath += "/"
if ScriptPath == "none/":
	ScriptPath = ""
OutFolder = sys.argv[11]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[12]
if OutFilePre == "none":
	OutFilePre = ""
PFileName = sys.argv[13]
NCores = sys.argv[14]
Mode = sys.argv[15]
if Mode not in ModeList:
	sys.exit("ERROR!  You wanted to use the mode %s, but it can only be one of the following: %s.\n%s" % (Mode, ", ".join(ModeList), Usage))
if Mode == "Array":
	NextScript = sys.argv[16]

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
	
#ListDictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value.  Each key has multiple values, and they
#will be made into a list
#modified from DictFromFile
def ListDictFromFile(FileName):
	TempDict = defaultdict(list)
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]].append(Line[1])
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is LPDict

#LocusSeqCombiner reads a series of sequence files and makes a dictionary of the sequences
#that have been classified according to locus.  Modified from the version in tbaits_intron_removal.py
#Modified from tseq_placer.py
def LocusSeqCombiner(LList,GDict,FilePath,FolderPre,FolderPost,FilePre,FilePost,SeqFormat):
	NumFiles = 0
	TempDict = defaultdict(dict)#TempDict[Paralog][SeqName] = ParSeq
	TempNameDict = { }
	for Locus in LList:
		ParName = Locus.split("_P")[0]
		for Group in GDict:
			FileName = FilePath+Group+"/"+FolderPre+GDict[Group]+"/"+FolderPost+FilePre+Locus+FilePost
			if (FolderPre == "") and (GDict[Group] == ""):
				FileName = FilePath+Group+"/"+FolderPost+FilePre+Locus+FilePost
			try:
				InFile = open(FileName, 'rU')
				for record in SeqIO.parse(InFile, SeqFormat):
					TempDict[ParName][record.id] = str(record.seq)
					TempNameDict[record.id] = ParName
				#print("File %s for locus %s was read.\n" % (FileName, Locus))
				InFile.close()
				NumFiles += 1
			except IOError:
				"alas, no file :("
				#print("File %s was not found.\n" % (FileName))
		#print("%d sequences were found for locus %s.\n" % (len(TempDict[Locus]), Locus))
	#print("%d sequence files for %d paralogs were read, with names such as %s.\n" % (NumFiles, len(LList), FileName))
	#sys.stderr.write("%d sequence files for %d paralogs were read, with names such as %s.\n" % (NumFiles, len(LList), FileName))
	return TempDict, TempNameDict
	#These are ParDict and LPNameDict

#LocusSeqListMaker reads a file of sequences and makes a list of the sequence names in that file
#original
def LocusSeqListMaker(FileName,SeqFormat):
	ListTemp = [ ]
	InFile = open(FileName, 'rU')
	for record in SeqIO.parse(InFile, SeqFormat):
		ListTemp.append(record.id)
	InFile.close()
	return ListTemp
	#This is LocusSeqList

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

#################################################################################


#make the dictionary of the taxonomic groups and their shortened names
#GroupDict[Group] = GroupAbb
GroupDict = NoneDictFromFile(GFileName)
#make the dictionary of outgroups
OutGroupDict = DictFromFile(OGFileName)

#reading the information files and condensing them:
#First: OutFolder+OutFilePre_Ind_Seq_Info.txt (tab-delimitted)
'''[0] Locus: alaAT
[1] Paralog:Ambig_alaAT2_alaAT1
[2] IndName: Alluaudia_procera_53
[3] #_Sequences: 4
[4] #_Ambiguities:3 
[5] #Seq_Length: 507
[6] Using_Seq: no'''
IndSeqInfoDict = defaultdict(dict)#IndSeqInfoDict[Locus][Paralog][IndName] = info
IndAmbigDict = defaultdict(dict)#IndAmbigDict[Locus][Paralog][IndName] = % ambiguities
GroupIndDict = defaultdict(list)
LPDict = defaultdict(list)
GroupList = GroupDict.keys()
for Group in GroupList:
	InFileName = SeqFilePath+Group+"/"+SeqFolderPre+GroupDict[Group]+"/"+SeqFolderPost+SeqFilePre+"Ind_Seq_Info.txt"
	if (SeqFolderPre == "") and (GroupDict[Group] == ""):
		InFileName = SeqFilePath+Group+"/"+SeqFolderPost+SeqFilePre+"Ind_Seq_Info.txt"
	try:
		InFile = open(InFileName, 'rU')
		for Line in InFile:
			Line = Line.strip('\r').strip('\n').split('\t')
			if Line[0] != "Locus":
				Locus = Line[0]
				Paralog = Line[1]
				IndName = Line[2]
				GroupIndDict[Group].append(IndName)
				try:
					IndSeqInfoDict[Locus][Paralog][IndName] = Line[3]+" seqs, overlap: "+Line[4]+", "+Line[5]+" ambigs, length: "+Line[6]+", using: "+Line[7]
					if (Line[5] != "n/a") and (Line[6] != "n/a"):
						IndAmbigDict[Locus][Paralog][IndName] = float(Line[5])/float(Line[6])
					else:
						IndAmbigDict[Locus][Paralog][IndName] = "n/a"
				except KeyError:
					IndSeqInfoDict[Locus]= defaultdict(dict)
					IndSeqInfoDict[Locus][Paralog][IndName] = Line[3]+" seqs, overlap: "+Line[4]+", "+Line[5]+" ambigs, length: "+Line[6]+", using: "+Line[7]
					IndAmbigDict[Locus] = defaultdict(dict)
					if (Line[5] != "n/a") and (Line[6] != "n/a"):
						IndAmbigDict[Locus][Paralog][IndName] = float(Line[5])/float(Line[6])
					else:
						IndAmbigDict[Locus][Paralog][IndName] = "n/a"
		InFile.close()
	except IOError:
		del GroupDict[Group]
		print("Group %s will not be examined further, as it does not appear to have been in this round of sequence analysis.\n" % (Group))
		sys.stderr.write("Group %s will not be examined further, as it does not appear to have been in this round of sequence analysis.\n" % (Group))
for Group in GroupIndDict:
	ListTemp = GroupIndDict[Group]
	GroupIndDict[Group] = sorted(list(set(ListTemp)))
for Locus in IndSeqInfoDict:
	for Paralog in IndSeqInfoDict[Locus]:
		LPDict[Locus].append(Paralog)

IndList = [ ]
for Group in GroupIndDict:
	IndList += GroupIndDict[Group]
OutFileName = OutFolder+OutFilePre+"Ind_Seq_Info.txt"
OutList = ["Locus\tParalog\t"+'\t'.join(IndList)]
for Locus in IndSeqInfoDict:
	for Paralog in IndSeqInfoDict[Locus]:
		Line = "\n"+Locus+"\t"+Paralog
		for IndName in IndList:
			try:
				Line += "\t"+IndSeqInfoDict[Locus][Paralog][IndName]
			except KeyError:
				Line += "\t"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)
OutFileName = OutFolder+OutFilePre+"perc_ambig.txt"
OutList = ["Locus\tParalog\t"+'\t'.join(IndList)]
for Locus in IndAmbigDict:
	for Paralog in IndAmbigDict[Locus]:
		Line = "\n"+Locus+"\t"+Paralog
		for IndName in IndList:
			try:
				if IndAmbigDict[Locus][Paralog][IndName] != "n/a":
					Line += ("\t%.4f" % (IndAmbigDict[Locus][Paralog][IndName]))
				else:
					Line += "\tn/a"
			except KeyError:
				Line += "\tn/a"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

OutFileName = OutFolder+OutFilePre+"Locus_Paralog_List.txt"
OutList = [ ]
for Locus in LPDict:
	for Paralog in LPDict[Locus]:
		Line = Locus+"\t"+Paralog+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#Second: OutFolder+OutFilePre_Paralog_Info.txt (tab-delimitted)
'''[0] Locus: ppc
[1] Paralog: ppc-1E2
[2] #_Ambiguities_in_Good_Seqs: 0,0,0,0,0,0,0,0,0 (comma-delimitted, can be blank if there are no good sequences)
[3] #_Ambiguities_in_Bad_Seqs: 13,21 (comma-delimitted, can be blank if there are no bad sequences)
[4] #_Good_Seqs: 9
[5] #_Bad_Seqs: 2'''
GroupInfoDict = defaultdict(dict)
for Group in GroupDict:
	InFileName = SeqFilePath+Group+"/"+SeqFolderPre+GroupDict[Group]+"/"+SeqFolderPost+SeqFilePre+"Paralog_Info.txt"
	if (SeqFolderPre == "") and (GroupDict[Group] == ""):
		InFileName = SeqFilePath+Group+"/"+SeqFolderPost+SeqFilePre+"Paralog_Info.txt"
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] != "Locus":
			Locus = Line[0]
			Paralog = Line[1]
			try:
				GroupInfoDict[Locus][Paralog][Group+"_Ambig"] = Line[2]+"\t"+Line[3]
			except KeyError:
				GroupInfoDict[Locus]= defaultdict(dict)
				GroupInfoDict[Locus][Paralog][Group+"_Ambig"] = Line[2]+"\t"+Line[3]
			GroupInfoDict[Locus][Paralog][Group+"_S"] = Line[4]+"\t"+Line[5]
	InFile.close()
OutFileName1 = OutFolder+OutFilePre+"Ambigs_per_Seq.txt"
OutFileName2 = OutFolder+OutFilePre+"Seqs_per_Locus.txt"
GroupList = sorted(GroupDict.keys())
OutString1 = "Locus\tParalog"
OutString2 = "Locus\tParalog"
for Group in GroupList:
	OutString1 += "\t"+Group+"_#_Ambigs_in_Good_Seqs\t"+Group+"_#_Ambigs_in_Bad_Seqs"
	OutString2 += "\t"+Group+"_Total_Good_Seqs\t"+Group+"_Total_Bad_Seqs"
OutList1 = [OutString1]
OutList2 = [OutString2]
for Locus in GroupInfoDict:
	for Paralog in GroupInfoDict[Locus]:
		Line1 = "\n"+Locus+"\t"+Paralog
		Line2 = "\n"+Locus+"\t"+Paralog
		for Group in GroupList:
			try:
				Line1 += "\t"+GroupInfoDict[Locus][Paralog][Group+"_Ambig"]
				Line2 += "\t"+GroupInfoDict[Locus][Paralog][Group+"_S"]
			except KeyError:
				Line1 += "\t\t"
				Line2 += "\t\t"
		OutList1.append(Line1)
		OutList2.append(Line2)
OutFileWriting(OutFileName1, OutList1)
OutFileWriting(OutFileName2, OutList2)

#Third, OutFolder+OutFilePre+Contig_Fates.txt
ContigFateDict = defaultdict(dict)
for Group in GroupDict:
	InFileName = SeqFilePath+Group+"/"+SeqFolderPre+GroupDict[Group]+"/"+SeqFolderPost+SeqFilePre+"Contig_Fates.txt"
	if (SeqFolderPre == "") and (GroupDict[Group] == ""):
		InFileName = SeqFilePath+Group+"/"+SeqFolderPost+SeqFilePre+"Contig_Fates.txt"
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] != "Locus":
			Locus = Line[0]
			Contig = Line[1]
			Fate = Line[2]
			try:
				ContigFateDict[Locus][Group][Contig] = Fate
			except KeyError:
				ContigFateDict[Locus][Group] = defaultdict(dict)
				ContigFateDict[Locus][Group][Contig] = Fate
	InFile.close()
OutFileName = OutFolder+OutFilePre+"Contig_Fates.txt"
OutList = ["Locus\tGroup\tContigs\tFate\n"]
for Locus in sorted(ContigFateDict.keys()):
	for Group in sorted(ContigFateDict[Locus].keys()):
		for Contig in sorted(ContigFateDict[Locus][Group].keys()):
			Line = Locus+"\t"+Group+"\t"+Contig+"\t"+ContigFateDict[Locus][Group][Contig]+"\n"
			OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#finding the loci and writing them to the files
LPNameDict = defaultdict(dict)
OutLocusDict = defaultdict(list) 
SeqsperLocus = { }
FileMovingScriptName = OutFolder+OutFilePre+"file_moving_script.sh"
FileMovingScript = ['#! /bin/bash\n\n']
for Locus in LPDict:
	(ParDict, LPNameDict[Locus]) = LocusSeqCombiner(LPDict[Locus], GroupDict, SeqFilePath, SeqFolderPre, SeqFolderPost, SeqFilePre, ".fa", "fasta")
	#LocusSeqCombiner(LList,GDict,FilePath,FolderPre,FolderPost,FilePre,FilePost,SeqFormat)
	#FileName = FilePath+Group+"/"+FolderPre+GroupDict[Group]+FolderPost+FilePre+Locus+FilePost
	#FileName = SeqFilePath+Group+"/"+SeqFolderPre+GroupDict[Group]+SeqFolderPost+SeqFilePre+Locus+".fa"
	LocusSeqList = LocusSeqListMaker(AlFolder+AlFilePre+Locus+AlFilePost+".fa", "fasta")
	if len(ParDict.keys()) != 0:
		NumSeqs = 0
		NumGoodSeqs = 0
		#making one file that contains all of the sequences for that locus
		OutFileName1 = OutFolder+OutFilePre+Locus+"_combined_all.fa"
		OutFile1 = open(OutFileName1, 'w')
		#making another file that contains only the good sequences for that locus
		OutFileName3 = OutFolder+OutFilePre+Locus+"_combined_best.fa"
		OutFile3 = open(OutFileName3, 'w')
		for Paralog in ParDict:
			SeqLenList = [ ]
			NewParDict = { }
			for ParName in ParDict[Paralog]:
				SeqLenList.append(len(ParDict[Paralog][ParName]))
			SeqMeanLen = numpy.mean(SeqLenList)
			#and making a third (bzw. 2nd) outfile that only contains the sequences from a single paralog
			OutFileName2 = OutFolder+OutFilePre+Locus+"_"+Paralog+".fa"
			OutFile2 = open(OutFileName2, 'w')
			for ParName in ParDict[Paralog]:
				if ParName in LocusSeqList:
					print("Sequence %s is already present in the alignment!\n" % (ParName))
					SplitParName = ParName.split(".")
					#Pereskia_saccharosa_67.10951.1seqs.0ambig.len818
					LocusPart = SplitParName[1]+"_new"
					SplitParName[1] = LocusPart
					NewParName = ".".join(SplitParName)
					if NewParName in LocusSeqList:
						LocusPart += "b"
						SplitParName[1] = LocusPart
						NewParName = ".".join(SplitParName)
						if NewParName in LocusSeqList:
							sys.exit("ERROR!!  The name %s has been repeated too many times; there is likely a problem!!\n" % (ParName))
					print("The new name for this sequence is %s.\n" % (NewParName))
					NewParDict[NewParName] = ParName
			#if any repeated names were found,
			if NewParDict != { }:
				for NewParName in NewParDict:
					ParName = NewParDict[NewParName]
					#add the new name to ParDict
					ParDict[Paralog][NewParName] = ParDict[Paralog][ParName]
					#and delete the old name
					print("Sequence %s will be removed from the ParDict.\n" % (ParName))
					del ParDict[Paralog][ParName]
			#then we can go back to writing sequences
			for ParName in ParDict[Paralog]:
				#print ("%s: %d." % (ParName, len(ParDict[Paralog][ParName])))
				Record1 = SeqRecord(seq=Seq(ParDict[Paralog][ParName], IUPAC), id = ParName, description = "")
				#definitely add the sequence to the file for that paralog, regardless of how long it is
				SeqIO.write(Record1, OutFile2, "fasta")
				#add it to the file containing all sequences for that locus regardless of how long it is
				SeqIO.write(Record1, OutFile1, "fasta")
				#if it is better, add it to the file containing the best sequences for that locus, the ones that will be used for the new backbone tree
				#*********We currently accept all sequences that are at least 3/4 the mean sequence length and at least 75 bp, but this can be changed.********
				if (len(ParDict[Paralog][ParName]) > 0.75*SeqMeanLen) and (len(ParDict[Paralog][ParName]) > 75):
					SeqIO.write(Record1, OutFile3, "fasta")
					NumGoodSeqs += 1
					#print ("using!")
				NumSeqs += 1
				#print("Locus: %s, Paralog: %s, Sequence Name: %s, Sequence Length: %d, Sequence Number: %d\n" % (Locus, Paralog, ParName, len(ParDict[Paralog][ParName]), NumSeqs))
			OutFile2.close()
		OutFile1.close()
		print("%d sequences were written to the file %s.\n" % (NumSeqs, OutFileName1))
		#sys.stderr.write("%d sequences were written to the file %s.\n" % (NumSeqs, OutFileName1))
		print("%d of these sequences were long enough to be written to the new backbone file %s.\n" % (NumGoodSeqs, OutFileName3))
		if NumGoodSeqs != 0:
			OutLocusDict[NumGoodSeqs].append(Locus)
			SeqsperLocus[Locus] = NumGoodSeqs
		else:
			Line = "cp "+AlFolder+AlFilePre+Locus+AlFilePost+".fa "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
			Line += "cp "+AlFolder+"RAxML_bestTree."+AlFilePre+Locus+" "+OutFolder+"RAxML_bestTree."+OutFilePre+Locus+"\n"
			FileMovingScript.append(Line)
		#sys.stderr.write("%d of these sequences were long enough to be written to the new backbone file %s.\n" % (NumGoodSeqs, OutFileName3))
		print("These sequences were also written to %d separate files for each paralog, with names such as %s.\n" % (len(ParDict), OutFileName2))
	else:
		Line = "cp "+AlFolder+AlFilePre+Locus+AlFilePost+".fa "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		Line += "cp "+AlFolder+"RAxML_bestTree."+AlFilePre+Locus+" "+OutFolder+"RAxML_bestTree."+OutFilePre+Locus+"\n"
		FileMovingScript.append(Line)
#writing the script to analyze the files for the bad loci
OutFileWriting(FileMovingScriptName, FileMovingScript)

#making additions to the pardict, the file that says which sequence belongs to which paralog
OutList = [ ]
#Then adding the new sequences to the pardict.  They need to be of the format:
#Locus [tab] SeqName [tab] Paralog [tab] Genus
for Locus in sorted(LPNameDict.keys()):
	for SeqName in LPNameDict[Locus]:
		Genus = SeqName.split("_")[0]
		Line = Locus+"\t"+SeqName+"\t"+LPNameDict[Locus][SeqName]+"\t"+Genus+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"pardict.txt"
OutFileWriting(OutFileName, OutList)


OutLocusList = [ ]
for NumGoodSeqs in sorted(OutLocusDict.keys(), reverse=True):
	if NumGoodSeqs != 0:
		OutLocusList += OutLocusDict[NumGoodSeqs]

if Mode == "Parallel":
	#the overall analysis script
	OutList1 = ["#! /bin/bash\n"]
	OutFileName1 = OutFolder+OutFilePre+"Analysis_Script.sh"
	#the script that can be parallelized
	OutList2 = ["#! /bin/bash\n"]
	OutFileName2 = OutFolder+OutFilePre+"Al_Tr_Script.sh"
	#first, add the old pardict to the new pardict; we want to do it this way, instead of just automatically writing the pardict
	#so that mutliple rounds of tparalog_combiner.py can be combined by another script more easily.
	Line = "cat "+OutFolder+OutFilePre+"pardict.txt "+PFileName+" > "+OutFolder+OutFilePre+"pardict_new.txt\n"
	#then change the mode of the second file, so we can execute it.
	Line += "chmod u+x "+OutFileName2+"\n"
	OutList1.append(Line)
	#then deal with the alignments
	for Locus in OutLocusList:
		Line = "rm "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		Line += "rm "+OutFolder+"RAxML*"+OutFilePre+Locus+"\n"
		OutList1.append(Line)
		Line = "mafft --localpair --add "+OutFolder+OutFilePre+Locus+"_combined_best.fa --quiet --thread -1 "+AlFolder+AlFilePre+Locus+AlFilePost+".fa > "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa && "
		Line += ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa && "
		#Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+"_allseqs_al.phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
		Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+"_allseqs_al.phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f d -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
		OutList2.append(Line)
	Line = "cat "+OutFileName2+" | parallel --jobs "+NCores+" --joblog "+OutFolder+OutFilePre+"parallel_log.log\n"
	OutList1.append(Line)
	OutFileWriting(OutFileName1, OutList1)
	OutFileWriting(OutFileName2, OutList2)
	print("The shell script for analyzing the sequences further was written to %s.\n" % (OutFileName1))
	sys.stderr.write("The shell script for analyzing the sequences further was written to %s.\n" % (OutFileName1))
elif Mode == "Array":
	SBatchList = [ ]
	#each gene will get its own file (although this is not actually an array)
	for Locus in OutLocusList:
		OutScript = ["#!/bin/bash\n"]
		OutScript.append("#SBATCH -J "+OutFilePre+Locus+"\n")
		OutScript.append("#SBATCH -t 6:00:00\n")
		OutScript.append("#SBATCH -n 1\n")
		#giving more cores to large alignments
		if (SeqsperLocus[Locus] > 50) or (Locus in ["ppc1", "ppc2"]):
			OutScript.append("#SBATCH --mem=16G\n")
		else:
			OutScript.append("#SBATCH --mem=4G\n")
		Line = "module load mafft\nmodule load raxml\n"
		OutScript.append(Line)
		Line = "rm "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		Line += "rm "+OutFolder+"RAxML*"+OutFilePre+Locus+"\n"
		OutScript.append(Line)
		Line = "mafft --localpair --add "+OutFolder+OutFilePre+Locus+"_combined_best.fa --quiet --thread -1 "+AlFolder+AlFilePre+Locus+AlFilePost+".fa > "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		Line += ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		#I don't think we use the bootstrap values, so I can save time by just making a tree.
		Line += "raxmlHPC -f d -s "+OutFolder+OutFilePre+Locus+"_allseqs_al.phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
		#Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+"_allseqs_al.phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
		OutScript.append(Line)
		OutFileName = OutFolder+OutFilePre+Locus+"_script.sh"
		OutFileWriting(OutFileName, OutScript)
		OutLine = "sbatch "+OutFileName
		SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		SBatchList.append(SBatchOut.strip('\r').strip('\n').split(" ")[-1])
	print("%d separate jobs were submitted to reconstruct the gene trees.\n" % (len(SBatchList)))
	sys.stderr.write("%d separate jobs were submitted to reconstruct the gene trees.\n" % (len(SBatchList)))
	OutScript = ['#! /bin/bash\n#SBATCH -J '+OutFilePre+'parcombiner\n#SBATCH -t 2:00:00\n\n']
	#first, add the old pardict to the new pardict; we want to do it this way, instead of just automatically writing the pardict
	#so that mutliple rounds of tparalog_combiner.py can be combined by another script more easily.
	Line = "cat "+OutFolder+OutFilePre+"pardict.txt "+PFileName+" > "+OutFolder+OutFilePre+"pardict_new.txt\n"
	OutScript.append(Line)
	#running the file moving script:
	Line = "chmod u+x "+FileMovingScriptName+"\n"
	Line += FileMovingScriptName+"\n"
	OutScript.append(Line)
	#running the next master script:
	Line = "chmod u+x "+NextScript+'\n'
	Line += NextScript+'\n'
	OutScript.append(Line)
	OutFileName = OutFolder+OutFilePre+"Analysis_Script.sh"
	OutFileWriting(OutFileName, OutScript)
	#this will take place only after the job arrays for all of the individuals have been run.
	if len(SBatchList) > 0:
		OutLine = "sbatch -d afterok:"+":".join(SBatchList)+" "+OutFileName
	#But it is possible that everything has been run, in which case it will just make the next files.
	else:
		OutLine = "sbatch "+OutFileName
	SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
	JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
	print("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
	sys.stderr.write("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
