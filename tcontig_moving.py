#! /usr/bin/env python

#tcontig_moving.py version 1.0 12 April 2016 Abby Moore
#This script moves contigs from the group folders where they were made to folders in groups that make more sense for tree building

import sys#getting information from the command line
import subprocess#We want to be able to get input from shell processes
from collections import defaultdict#dictionaries with multiple levels
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

'''
tcontig_moving.py InFilePath OutFilePath OldGroupPreFN NewGroupPreFN IndGroupFN LocusListFN Mode[Spades, Minimo, Seqs]
/users/ajm3/data/ajm3/scripts/tcontig_moving.py /gpfs/scratch/ajm3/eedwards/ /gpfs/scratch/ajm3/eedwards/ /users/ajm3/data/ajm3/general/GroupList_Ln12345.txt /gpfs/scratch/ajm3/eedwards/comb_p1_spades/Comb_paper1_Group_List.txt /gpfs/scratch/ajm3/eedwards/comb_p1_spades/comb_Anac.txt /gpfs/scratch/ajm3/eedwards/general/Locus_List_alaAT.txt Spades
'''

Usage = '''
tcontig_moving.py
This script moves contigs from the group folders where they were made to folders
that make more sense for sequence making and tree building.
[path to the groups of input files]
[path to where the output files should go]
[file with OldGroupName [tab] OldGroupPrefix]
[file with NewGroupName [tab] NewGroupPrefix]
[file with OldGroup [tab] IndName [tab] NewGroup]
[file with the list of loci]
[Mode: Spades, Minimo, or Seqs]
'''

ModeList = ['Spades', 'Minimo', 'Seqs']

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 8:
	sys.exit("ERROR!!  This script requires at least 7 additional arguments and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFilePath = sys.argv[1]
if InFilePath[-1] != "/":
	InFilePath += "/"
OutFilePath = sys.argv[2]
if OutFilePath[-1] != "/":
	OutFilePath += "/"
OldGroupPreFN = sys.argv[3]
NewGroupPreFN = sys.argv[4]
IndGroupFN = sys.argv[5]
LocusListFN = sys.argv[6]
Mode = sys.argv[7]
if Mode not in ModeList:
	sys.exit("ERROR!!  You wanted mode %s, but it must be one of the following: %s.\n%s" % (Mode, ", ".join(ModeList), Usage))


Vociferous = False

Vociferous = True

#################################################################################

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
	#This is various things

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

#SeqFileReading reads a sequence file and puts the sequences in a dictionary.
#from tparcomb_final.py (although it was maybe from somewhere else before that)
def SeqFileReading(FileName, SeqFormat):
	DictTemp = { }#DictTemp[SeqName] = Seq
	InFile = open(FileName, 'rU')
	for record in SeqIO.parse(InFile, SeqFormat):
		DictTemp[record.id] = str(record.seq)
	InFile.close()
	return DictTemp

#SeqFileWriting writes sequence files from dictionaries.
#from tparcomb_final.py (although it was maybe from somewhere else before that)
def SeqFileWriting(FileName, SDict, SeqFormat):
	OutFile = open(FileName, 'w')
	for SeqName in SDict:
		Record1 = SeqRecord(seq=Seq(SDict[SeqName], IUPAC), id = SeqName, description = "")
		SeqIO.write(Record1, OutFile, SeqFormat)
	OutFile.close()


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

#making dictionaries from the input files
OldGroupDict = DictFromFile(OldGroupPreFN, 0, 1)#OldGroupDict[OldGroupName] = OldGroupPre
NewGroupDict = DictFromFile(NewGroupPreFN, 0, 1)#NewGroupDict[NewGroupName] = NewGroupPre
#IndOldGroups = DictFromFile(IndGroupFN, 0, 1)#IndOldGroups[IndName] = OldGroupName
#IndNewGroups = DictFromFile(IndGroupFN, 0, 2)#IndNewGroups[IndName] = NewGroupName
OGIndList = ListDictFromFile(IndGroupFN, 0, 1)#OGIndList[OldGroupName] = [list of IndNames]
NGIndList = ListDictFromFile(IndGroupFN, 2, 1)#NGIndList[NewGroupName] = [list of IndNames]
LocusList = CaptureColumn(LocusListFN, 0)#list of loci

if Mode == "Spades":
	IntFolderPre = "s2_final_"
	FinalFolder = "spades_contigs"
	ContigFilePre = "sc_"
	BlastFilePre = "sb3_"
elif Mode == "Minimo":
	IntFolderPre = "min_final_"
	FinalFolder= "minimo_contigs"
	ContigFilePre = "minc_"
	BlastFilePre = "minb3_"
elif Mode == "Seqs":
	IntFolderPre = "sts_final_"
	FinalFolder = "spades_contigs"
	SeqFilePre = "stsf_"

#first, make the new folders
#OutFileName = OutFilePath+Mode+"_file_moving.sh"
#OutList = ["#! /bin/bash\n#SBATCH -J "+Mode+"_file_moving\n#SBATCH -t 2:00:00\n\n\n"]
for NewGroupName in NewGroupDict:
	OutCode = subprocess.call("mkdir "+OutFilePath+NewGroupName, shell = True)
	OutCode = subprocess.call("mkdir "+OutFilePath+NewGroupName+"/"+IntFolderPre+NewGroupDict[NewGroupName], shell = True)
	OutCode = subprocess.call("mkdir "+OutFilePath+NewGroupName+"/"+IntFolderPre+NewGroupDict[NewGroupName]+"/"+FinalFolder, shell = True)
#then, make the individual contig files
#going locus by locus
if (Mode == "Spades") or (Mode == "Minimo"):
	for Locus in LocusList:
		NumSeqFiles = 0
		NumBlastFiles = 0
		#make a dictionary of all of the sequences from the individuals that we want from that locus
		LocusDict = defaultdict(dict)#LocusDict[IndName][SeqName] = Seq
		#do this by reading the old sequence files one by one
		for OldGroupName in OGIndList:
			InFileName = InFilePath+OldGroupName+"/"+IntFolderPre+OldGroupDict[OldGroupName]+"/"+FinalFolder+"/"+ContigFilePre+Locus+".fa"
			try:
				LocusDictTemp = SeqFileReading(InFileName, "fasta")
				#looking at each sequence name
				for SeqName in LocusDictTemp:
					IndName = SeqName.split("-")[0]
					#and adding it to the LocusDict, if that individual is one we want
					if IndName in OGIndList[OldGroupName]:
						LocusDict[IndName][SeqName] = LocusDictTemp[SeqName]
			except IOError:
				"do nothing"
		#then write them to the files, by going through each new group name
		for NewGroupName in NGIndList:
			LocusDictTemp = { }
			#and pulling out all of the sequences for the individuals in that group from the LocusDict
			for IndName in NGIndList[NewGroupName]:
				if IndName in LocusDict.keys():
					for SeqName in LocusDict[IndName]:
						LocusDictTemp[SeqName] = LocusDict[IndName][SeqName]
				else:
					if Vociferous == True: print("No sequences for individual %s for locus %s.\n" % (IndName, Locus))
			#and writing them to their new files
			OutFileName = OutFilePath+NewGroupName+"/"+IntFolderPre+NewGroupDict[NewGroupName]+"/"+FinalFolder+"/"+ContigFilePre+Locus+".fa"
			SeqFileWriting(OutFileName, LocusDictTemp, 'fasta')
			NumSeqFiles += 1
		#doing the same thing for the blast output
		LocusBlastDict = defaultdict(list)#LocusDict[IndName] = [list of Lines]
		for OldGroupName in OGIndList:
			InFileName = InFilePath+OldGroupName+"/"+IntFolderPre+OldGroupDict[OldGroupName]+"/"+FinalFolder+"/"+BlastFilePre+Locus+".out"
			try:
				InFile = open(InFileName, 'rU')
				for Line in InFile:
					Line = Line.strip('\n').strip('\r')
					#getting the individual name from the line
					SplitLine = Line.split('\t')
					IndName = SplitLine[0].split("-")[0]
					if IndName in OGIndList[OldGroupName]:
						LocusBlastDict[IndName].append(Line)
				InFile.close()
			except IOError:
				"do nothing"
		for NewGroupName in NGIndList:
			LocusBlastDictTemp = { }
			#and pulling out all of the sequences for the individuals in that group from the LocusBlastDict
			for IndName in NGIndList[NewGroupName]:
				if IndName in LocusBlastDict.keys():
					LocusBlastDictTemp[IndName] = LocusBlastDict[IndName]
				else:
					if Vociferous == True: print("No blast hits for individual %s for locus %s.\n" % (IndName, Locus))
			#and writing them to their new files
			OutFileName = OutFilePath+NewGroupName+"/"+IntFolderPre+NewGroupDict[NewGroupName]+"/"+FinalFolder+"/"+BlastFilePre+Locus+".out"
			OutFile = open(OutFileName, 'w')
			for IndName in LocusBlastDictTemp:
				for Line in LocusBlastDictTemp[IndName]:
					OutFile.write(Line+"\n")
			OutFile.close()
			NumBlastFiles += 1
		if NumSeqFiles != NumBlastFiles:
			sys.exit("ERROR!!!  The number of sequence and blast files should be the same, but there were %d sequence files and %d blast files!" % (NumSeqFiles, NumBlastFiles))
		print("%d sequence files and %d blast files were written for locus %s.\n" % (NumSeqFiles, NumBlastFiles, Locus))
		sys.stderr.write("%d sequence files and %d blast files were written for locus %s.\n" % (NumSeqFiles, NumBlastFiles, Locus))
elif Mode == "Seqs":
	for Locus in LocusList:
		NumSeqFiles = 0
		#make a dictionary of all of the sequences from the individuals that we want from that locus
		LocusDict = defaultdict(dict)#LocusDict[IndName][SeqName] = Seq
		#do this by reading the old sequence files one by one
		for OldGroupName in OGIndList:
			InFileName = InFilePath+OldGroupName+"/"+IntFolderPre+OldGroupDict[OldGroupName]+"/"+FinalFolder+"/"+SeqFilePre+Locus+".fa"
			try:
				LocusDictTemp = SeqFileReading(InFileName, "fasta")
				#looking at each sequence name
				for SeqName in LocusDictTemp:
					IndName = SeqName.split("-")[0]
					#and adding it to the LocusDict, if that individual is one we want
					if IndName in OGIndList[OldGroupName]:
						LocusDict[IndName][SeqName] = LocusDictTemp[SeqName]
						#if Vociferous == True: print("Sequence %s from individual %s in group %s will be used.\n" % (SeqName, IndName, OldGroupName))
					else:
						if Vociferous == True: print("Sequence %s from individual %s in group %s will NOT BE used.\n" % (SeqName, IndName, OldGroupName))
			except IOError:
				"do nothing"
		#then write them to the files, by going through each new group name
		for NewGroupName in NGIndList:
			LocusDictTemp = { }
			#and pulling out all of the sequences for the individuals in that group from the LocusDict
			for IndName in NGIndList[NewGroupName]:
				if IndName in LocusDict.keys():
					for SeqName in LocusDict[IndName]:
						LocusDictTemp[SeqName] = LocusDict[IndName][SeqName]
				else:
					if Vociferous == True: print("No sequences for individual %s for locus %s.\n" % (IndName, Locus))
			#and writing them to their new files
			OutFileName = OutFilePath+NewGroupName+"/"+IntFolderPre+NewGroupDict[NewGroupName]+"/"+FinalFolder+"/"+SeqFilePre+Locus+".fa"
			SeqFileWriting(OutFileName, LocusDictTemp, 'fasta')
			NumSeqFiles += 1
