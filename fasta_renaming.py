#!/usr/bin/env python

#fasta_renaming.py version 1.0 31 August 2015 Abby Moore

#This script renames sequences in fasta format (aligned or unaligned) according to a dictionary.
#The names must start with a certain phrase and then they will be given a new prefix.
#Names that do not have a prefix in the dictionary are kept but printed to the screen.

Usage = '''
fasta_renaming.py renames sequences in fasta format
[input folder]
[output folder, or same, if the same as the input folder]
[file of input files]
[new prefix for output files, or none, if the same name as the input files]
[tab-delimitted file with old prefix tab new prefix]
[format of the input sequences]
[desired format for the output sequences--However, it is not writing an alignment,
so it will not format nexus or phylip properly--best to make this in fasta and then
convert it]
'''

#fasta_renaming.py InFolder OutFolder ListFileName OutFilePre PreDictFileName SeqFormatIn SeqFormatOut

import sys
from collections import defaultdict
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Alphabet import IUPAC #to recognize sequences

if len(sys.argv) != 8:
	sys.exit("ERROR!  This script requires 7 additional arguments and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
if InFolder[-1] != "/":
	InFolder += "/"
OutFolder = sys.argv[2]
if OutFolder == "same":
	OutFolder = InFolder
elif OutFolder[-1] != "/":
	OutFolder += "/"
ListFileName = sys.argv[3]
OutFilePre = sys.argv[4]
if OutFilePre == "none":
	OutFilePre = ""
PreDictFileName = sys.argv[5]
SeqFormatIn = sys.argv[6]
SeqFormatOut = sys.argv[7]

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
#from tseq_placer.py
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
	#This is InFileList

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
	#This is PreDict.

#LocusSeqGetter reads a series of sequence files and makes a dictionary of the sequences
#that have been classified according to locus.  Used first in tbaits_intron_removal.py
#modified from tcontigs_to_fixed_paralogs.py
def LocusSeqGetter(FileList,Folder,SeqFormat):
	TempDict = defaultdict(dict)
	for SeqFile in FileList:
		FileName = Folder + SeqFile
		InFile = open(FileName, 'rU')
		for record in SeqIO.parse(InFile, SeqFormat):
			TempDict[SeqFile][record.id] = str(record.seq)
		print("%d sequences for locus %s were read from the file %s.\n" % (len(TempDict[SeqFile].keys()), SeqFile, FileName))
	print("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	sys.stderr.write("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	return TempDict
	#This is SeqDictIn

#SeqRenaming renames sequences according to a dictionary of old prefixes.
def SeqRenaming(SeqDict, NameDict):
	TempDict = defaultdict(dict)
	for SeqFile in SeqDict:
		print SeqFile
		for SeqName in SeqDict[SeqFile]:
			SeqTemp = SeqDict[SeqFile][SeqName]
			for Pre in NameDict:
				SeqFound = False
				if SeqName[:len(Pre)] == Pre:
					SeqFound = True
					NameTemp = NameDict[Pre]+SeqName
					TempDict[SeqFile][NameTemp] = SeqTemp
					#print "Found!!"
					break
			if SeqFound == False:
				TempDict[SeqFile][SeqName] = SeqTemp
				#print "Not Found!!
				print SeqName
	return TempDict
	#This is SeqDictOut

#LocusSeqWriter, modified from LocusParalogSeqWriter
#modified from tcontigs_to_fixed_paralogs.py
def LocusSeqWriter(SeqDict, Folder, Prefix, Suffix, SeqFormat):
	#SeqDict[Locus][ContigName] = Sequence
	TempDict = { }#TempDict[Locus] = file name
	for Locus in SeqDict:
		OutFileName = Folder+Prefix+Locus+Suffix
		TempDict[Locus] = Prefix+Locus
		OutFile = open(OutFileName, 'w')
		for ContigName in SeqDict[Locus]:
			Record1 = SeqRecord(seq=Seq(SeqDict[Locus][ContigName], IUPAC), id = ContigName, description = "")
			SeqIO.write(Record1, OutFile, SeqFormat)
		OutFile.close()
	print("%d files of renamed sequences were written, with names such as %s.\n" % (len(TempDict), OutFileName))
	sys.stderr.write("%d files of renamed sequences were written, with names such as %s.\n" % (len(TempDict), OutFileName))
	return TempDict


InFileList = CaptureColumn(ListFileName, 0)
PreDict = DictFromFile(PreDictFileName)
SeqDictIn = LocusSeqGetter(InFileList, InFolder, SeqFormatIn)
SeqDictOut = SeqRenaming(SeqDictIn, PreDict)
SeqNameDict = LocusSeqWriter(SeqDictOut, OutFolder, OutFilePre, "", SeqFormatOut)
