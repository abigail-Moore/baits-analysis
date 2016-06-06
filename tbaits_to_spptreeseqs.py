#! /usr/bin/env python

#tbaits_to_spptreeseqs.py version 1.0 17 April 2016 Abby Moore
#This script is based on tbaits_intron_removal.py but was rewritten for the backbone
#sequences for the original species tree (matK, rbcL, ndhF, ITS).  It determines 
#which parts of the sequences blast against the original sequences, but it does
#not make subsequent alignments, like tbaits_intron_removal.py does.
#LocusKeyFile:
'''PLocus [0] (apl1)
GLocus [1] (apl)'''


#examples:
'''

tbaits_to_spptreeseqs.py LocusKeyFile BlastFilePre SeqFilePre SeqFolder BlastFolder OutFolder OutFilePre 
/users/ajm3/data/ajm3/scripts/tbaits_to_spptreeseqs.py /users/ajm3/data/ajm3/general/Locus_List_spptree.txt stsb_ stsc_ /gpfs/scratch/ajm3/eedwards/Portulacaceae/sts_final_Port/spades_contigs/ same same stsf_
'''

from collections import defaultdict
import sys
from numpy import mean
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
import re#we need regular expressions to transfor the minimo sequences

Usage = '''
tbaits_to_spptreeseqs.py version 1.0
This script removes the ends sequences that don't blast to the correct species
tree backbone sequences.
tbaits_to_spptreeseqs.py
[tab-delimited file with paralog name (tab) gene family name]
[prefix for blast files, or "none", if none]
[prefix for sequence files, or "none", if none]
[folder where sequence files are found]
[folder where blast output files are found, or "same", if it is the same folder 
as where the sequences are found]
[folder where output files should be written, or "same", if it is the same 
folder as where the sequences are found]
[prefix for output files, or "none", if none]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 8:
	sys.exit("ERROR!  This script requires 7 additional arguments, and you supplied %d.  %s\n" % (len(sys.argv)-1, Usage))
LocusKeyFile = sys.argv[1]
BlastFilePre = sys.argv[2]
if BlastFilePre == "none":
	BlastFilePre = ""
SeqFilePre = sys.argv[3]
if SeqFilePre == "none":
	SeqFilePre = ""
SeqFolder = sys.argv[4]
BlastFolder = sys.argv[5]
if BlastFolder == "same":
	BlastFolder = SeqFolder
OutFolder = sys.argv[6]
if OutFolder == "same":
	OutFolder = SeqFolder
OutFilePre = sys.argv[7]
if OutFilePre == "none":
	OutFilePre = ""

#The lists and dictionaries we will fill out:
#BlastFileDict[BlastFile] = SeqFile (via DictFromFile)
#LocusDict[PLocus] = GLocus (via DictFromFile)
#SeqFileList list of all of the sequence file names (without their folder) (via ListFromDictValues)
#ContigDict[Locus][Contig] = ContigSeq (LocusSeqGetter)
#SeqPosDict[Locus][Contig][SeqPos] = #times that sequence was a found by blast (via SeqRangeDictMaker and BlastHitsCounter)
#ExonPosDict[Locus][Contig] = list of tuples with start and end of exons (via ExonFinder) 

MinRename = True

#adding slashes, if necessary:
if SeqFolder[-1] != "/":
	SeqFolder += "/"
if BlastFolder[-1] != "/":
	BlastFolder += "/"
if OutFolder[-1] != "/":
	OutFolder += "/"

###########read random files that give us information about the files we want to read

#DictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
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
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	print("Output file %s written.\n" % (FileName))
	sys.stderr.write("Output file %s written.\n" % (FileName))

#ListFromDictValues makes a shorted list of the values in a dictionary.
def ListFromDictValues(DictName):
	TempDict = { }
	for DictKey in DictName:
		TempDict[DictName[DictKey]] = DictKey
	TempList = sorted(TempDict.keys())
	return TempList

#SeqLenFinder finds the length of a sequence in an alignment, while ignoring all gaps
#from tnotung_homolog_parsing
def SeqLenFinder(Seq):
	Len = 0
	for Base in Seq:
		if Base != "-":
			Len += 1
	return Len

#LocusSeqGetter reads a series of sequence files and makes a dictionary of the sequences
#that have been classified according to locus.
def LocusSeqGetter(FileList,Folder,FilePre,FilePost,SeqFormat):
	TempDict = defaultdict(dict)
	for SeqFile in FileList:
		Locus = SeqFile.replace(FilePre,"").replace(FilePost,"")
		FileName = Folder + SeqFile
		try:
			InFile = open(FileName, 'rU')
			for record in SeqIO.parse(InFile, SeqFormat):
				TempDict[Locus][record.id] = str(record.seq)
				#print("%d sequences for locus %s were read from the file %s.\n" % (len(TempDict[Locus].keys()), Locus, FileName))
		except IOError:
			"alas, no file"
			print("The file %s was not found.\n" % (FileName))
			#sys.stderr.write("ERROR!!!\nERROR!!\nERROR!!\nThe file %s was not found!!\n" % (FileName))
	print("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	sys.stderr.write("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	return TempDict

#SeqRangeDictMaker makes a dictionary where each position in each contig (or other
#sequence has an entry, which is currently 0, but will become something else in the 
#next step).  It expects a dictionary of sequences classified by locus and sequenc name.
def SeqRangeDictMaker(DictName):
	TempDict = defaultdict(dict)
	for Locus in DictName:
		for SeqName in DictName[Locus]:
			SeqLen = len(DictName[Locus][SeqName])
			#print("%s: %s: %d\n" % (Locus,SeqName,SeqLen))
			TempDict[Locus][SeqName] = defaultdict(int)
			for SeqPos in range(1,SeqLen+1):
				TempDict[Locus][SeqName][SeqPos] = 0
			TempDict[Locus][SeqName]['R'] = 0
			TempDict[Locus][SeqName]['F'] = 0
	return TempDict

#BlastHitsCounter reads a set of blast files and counts how often the various bases
#are present in the contigs.  It expects the blast hits to be labeled according to
#paralog, with the paralog separated from the rest of the locus name by '+'.
#It then looks up the locus name in a dictionary.
def BlastHitsCounter(SeqDict, BlastList, Folder, eVal, LDict):
	for BlastFile in BlastList:
		try:
			InFile = open(Folder+BlastFile, 'rU')
			for Line in InFile:
				Line = Line.strip('\r').strip('\n').split('\t')
				eValue = float(Line[10])
				#if the blast hit is above the threshold e-value
				if eValue < eVal:
					Contig = Line[0]
					if MinRename == True:
						MinRe = re.compile("=")
						Contig = MinRe.sub("-",Contig) 
					BlastHit = Line[1].split('+')
					Locus = LDict[BlastHit[0]]
					if (int(Line[8]) > int(Line[9])):
						Direc = 'R'
					else:
						Direc = 'F'
					SeqRange = range(int(Line[6]),int(Line[7])+1)
					#print("%d-%d: %d: %s" % (int(Line[8]), int(Line[9]), len(SeqRange), Direc))
					for SeqPos in SeqRange:
						SeqDict[Locus][Contig][SeqPos] += 1
					SeqDict[Locus][Contig][Direc] += 1
		except IOError:
			"alas, no file"
			print("The file %s was not found.\n" % (Folder+BlastFile))
			#sys.stderr.write("The file %s was not found.\n" % (Folder+BlastFile))
	print("Blast files were read and exons with e-values greater than %e were retained.\n" % (eVal))
	sys.stderr.write("Blast files were read and exons with e-values greater than %e were retained.\n" % (eVal))
	return SeqDict

#ExonFinder looks at the dictionary produced by BlastHitsCounter and finds out which
#parts of the sequence have enough blast hits to count as exons.
def ExonFinder(SeqDict, MinNumHits):
	TempDict = defaultdict(dict)
	NumExons = 0
	for Locus in SeqDict:
		for Contig in SeqDict[Locus]:
			TempDict[Locus][Contig] = [ ]
			InExon = False
			for SeqPos in sorted(SeqDict[Locus][Contig].keys())[:-2]:
				#if SeqPos in ['F', 'R'] == False: #or can I just say all but the last two--I assume they will always be last
				if SeqDict[Locus][Contig][SeqPos] >= MinNumHits: #in the exon
					if InExon == False:
						#found the start of an exon
						ExonStart = SeqPos
						InExon = True
					#elif InExon == True: in the middle of the exon, do nothing
				else: #not in the exon
					if InExon == True: #found the end of an exon
						ExonEnd = SeqPos
						ExonTuple = (ExonStart, ExonEnd)
						TempDict[Locus][Contig].append(ExonTuple)
						InExon = False
						NumExons += 1
					#elif InExon == False: in the middle of an intron, do nothing
			if InExon == True: #if the contig is still in an exon when it finished
				ExonEnd = SeqPos
				ExonTuple = (ExonStart, ExonEnd)
				TempDict[Locus][Contig].append(ExonTuple)
				NumExons += 1
	print ("%d good sequences, eaach of which had at least %d blast hits, were found in the various contigs.\n" % (NumExons, MinNumHits))
	sys.stderr.write("%d good sequences, each of which had at least %d blast hits, were found in the various contigs.\n" % (NumExons, MinNumHits))
	return TempDict

#BlastHitSeqFinder takes a dictionary of sequences grouped by locus (output by LocusSeqGetter),
#the information about whether a sequence is in the same direction as the baits or 
#reversed (output by BlastHitsCounter), and information on exon position (output by
#ExonFinder) and makes new sequences with each intron replaced by two dashes.
#A separate dictionary is made for sequences that are sometimes match in a forward
#direction and sometimes match in a reverse direction.  These two dictionaries are
#combined into one for easier printing.
#modified from IntronPruner from tbaits_intron_removal.py
def BlastHitSeqFinder(SeqDict, FRDict, PosDict):
	TempDict = defaultdict(dict)
	NumSeqs = 0
	AmbigSeqs = 0
	NoExons = 0
	for Locus in SeqDict:
		for Contig in SeqDict[Locus]:
			IndName = Contig.split("-")[0]
			OldSeq = SeqDict[Locus][Contig]
			NewSeq = ""
			#print ("%s: %s" % (Locus,Contig))
			for ExonTuple in PosDict[Locus][Contig]:
				EStart = ExonTuple[0]-1
				EEnd = ExonTuple[1]
				NewSeq += OldSeq[EStart:EEnd]+"--"
			if NewSeq != "":
				NewSeq = NewSeq[:-2]
				NewName = Contig.split("-")[0]+"-"+Locus+"_"+Contig.split("-NODE_")[1][0]
				#if the sequence blasts in both directions ways, we don't want it.
				if (FRDict[Locus][Contig]['R'] != 0) and (FRDict[Locus][Contig]['F'] != 0):
					AmbigSeqs += 1
				#otherwise, we might want the sequence
				else: 
					if FRDict[Locus][Contig]['F'] == 0: #the contig is definitely reversed
						NewSeq2 = Seq(NewSeq, IUPAC.unambiguous_dna)
						NewSeq = str(NewSeq2.reverse_complement())
					NumSeqs += 1				
					try:
						TempDict[Locus][IndName][NewName] = NewSeq
					except KeyError:
						TempDict[Locus][IndName] = defaultdict(dict)
						TempDict[Locus][IndName][NewName] = NewSeq
			else:
				NoExons += 1
	print("%d unambiguous contigs and %d ambiguous contigs were transcribed.\n" % (NumSeqs, AmbigSeqs))
	sys.stderr.write("%d unambiguous contigs and %d ambiguous contigs were transcribed.\n" % (NumSeqs, AmbigSeqs))
	print("%d contigs did not blast strongly to the exon database.\n" % (NoExons))
	sys.stderr.write("%d contigs did not blast strongly to the exon database.\n" % (NoExons))
	return TempDict
	#TempDict is OutDict.

#LongestSeqFinder finds the longest sequence for each individual in a dictionary 
#of sequences of the format DictIn[Locus][IndName][ContigName] = Seq
#original
def LongestSeqFinder(SeqDictIn):
	TempDict = {}
	for Locus in SeqDictIn:
		TempDict[Locus] = defaultdict(dict)
		for IndName in SeqDictIn[Locus]:
			LongestSeqName = ""
			LongestSeqLen = 0
			for Contig in SeqDictIn[Locus][IndName]:
				ContigLen = len(SeqDictIn[Locus][IndName][Contig]) 
				if ContigLen > LongestSeqLen:
					LongestSeqLen = ContigLen
					LongestSeqName = Contig
			LongestSeq = SeqDictIn[Locus][IndName][LongestSeqName]
			SeqNameTemp = "_".join(LongestSeqName.split("_")[:-1])
			TempDict[Locus][SeqNameTemp] = LongestSeq
	return TempDict
	#This is OutDictLongest

#LocusSeqWriter writes sequences in the desired format to files from a
#dictionary where the top-level key is the group, the second-level key is the locus,
#the third-level key is the sequence name, and the value is the sequence.
#modified from LocusGroupSeqWriter from tbaits_intron_removal.py
def LocusSeqWriter(SeqDict, Folder, Prefix, Suffix, SeqSuffix, SeqFormat):
	TempDict = defaultdict(dict)
	NumFiles = 0
	for Locus in SeqDict:
		OutFileName = Folder+Prefix+Locus+Suffix
		TempDict[Locus] = Prefix+Locus
		OutFile = open(OutFileName, 'w')
		for Contig in SeqDict[Locus]:
			Record1 = SeqRecord(seq=Seq(SeqDict[Locus][Contig], IUPAC), id = Contig+SeqSuffix, description = "")
			SeqIO.write(Record1, OutFile, SeqFormat)
		OutFile.close()
		NumFiles += 1
	print("%d sequence files were written, with names such as %s.\n" % (NumFiles, OutFileName))
	sys.stderr.write("%d sequence files were written, with names such as %s.\n" % (NumFiles, OutFileName))
	return TempDict

#SeqStatsWriter reads the sequences and calculates some statistics.
#tbaits_intron_removal.py?
def SeqStatsWriter(SeqDict, Folder, Pre, Delimitter):
	TempDict = defaultdict(dict) #TempDict[Ind][Locus] = NumContigs
	MeanTempDict = { } #MeanTempDict[Locus][Ind] = [list of contig lengths]
	#making the dictionary showing the number of contigs per locus
	for Locus in SeqDict:
		MeanTempDict[Locus] = defaultdict(list)
		for ContigName in SeqDict[Locus]:
			IndName = ContigName.split(Delimitter)[0]
			MeanTempDict[Locus][IndName].append(SeqLenFinder(SeqDict[Locus][ContigName]))
			try:
				TempDict[IndName][Locus] += 1
			except KeyError:
				TempDict[IndName][Locus] = 1
	#writing the dictionary to the file (TempDict)
	TempIndList = sorted(TempDict.keys())
	OutList = ['Locus\t'+"\t".join(TempIndList)]
	for Locus in SeqDict:
		Line = "\n"+Locus
		for IndName in TempIndList:
			try:
				Line += "\t"+str(TempDict[IndName][Locus])
			except KeyError:
				Line += "\t"
		OutList.append(Line)
	OutFileName = Folder+Pre+"Contigs_per_Ind_Locus.txt"
	OutFileWriting(OutFileName, OutList)
	#writing MeanTempDict to a file
	OutList = ['Locus\t'+"\t".join(TempIndList)]
	for Locus in SeqDict:
		Line = "\n"+Locus
		for IndName in TempIndList:
			if MeanTempDict[Locus][IndName] != []:
				Line += "\t"+str(mean(MeanTempDict[Locus][IndName]))
			else:
				Line += "\t"
		OutList.append(Line)
	OutFileName = Folder+Pre+"Mean_Contig_Length_per_Ind_Locus.txt"
	OutFileWriting(OutFileName, OutList)
	return

#################################################################################

#Read the BlastFileList and the LocusKeyFile and make dictionaries from them:
LocusDict = DictFromFile(LocusKeyFile)
LocusList = ListFromDictValues(LocusDict)
SeqFileList = [SeqFilePre+Locus+".fa" for Locus in LocusList]
BlastFileList = [BlastFilePre+Locus+".out" for Locus in LocusList]

#read contig files, saving the contigs, and preparing the dictionary that will tell how many times blast found a given part of a sequence

ContigDict = LocusSeqGetter(SeqFileList, SeqFolder, SeqFilePre, ".fa", 'fasta')
SeqPosDicta = SeqRangeDictMaker(ContigDict)

#read blast files and figure out how many times each base in the contig file was hit

SeqPosDict = BlastHitsCounter(SeqPosDicta, BlastFileList, BlastFolder, 1e-8, LocusDict)

#determine which positions had enough hits to be exons

ExonPosDict = ExonFinder(SeqPosDict, 5)

#remove the exons and decide what to do with them when we reach that point

OutDict = BlastHitSeqFinder(ContigDict, SeqPosDict, ExonPosDict)

OutDictLongest = LongestSeqFinder(OutDict)

#Figure out which loci had sequences
SeqStatsWriter(OutDict, OutFolder, OutFilePre, "-")

#Write sequence files and shell scripts

FileDict = LocusSeqWriter(OutDictLongest, OutFolder, OutFilePre, '.fa', "", 'fasta')

