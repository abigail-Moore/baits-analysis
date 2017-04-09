#! /usr/bin/env python

#tbaits_intron_removal.py version 1.1 17 Nov. 2015 Abby Moore
#This script looks at the results from blasting genomic contigs against bait
#sequences and removes the intron sequences.  It is hopefully designed to be able
#to work with contigs made from gene families, where the exons might be in
#different places, unlike tbaits_to_exons.py
#version 1.1 has been changed to take LocusList as input and not rely on the BlastFileList
#File formats:
#BlastFileList:
'''BlastFileName [0] (bf2_nhd.out)
SeqFileName [1] (sf2_nhd.fa)'''
#LocusKeyFile:
'''PLocus [0] (apl1)
GLocus [1] (apl)'''


#examples:
'''
tbaits_intron_removal.py /home/abby/transcriptomes/TS31_1/sandbox/gfams/loci_shortened_ppc.txt bf2_ sf2_ ~/transcriptomes/TS31_1/sandbox/gfams/gfcontigs/ same ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/ of3t_ ~/transcriptomes/baits_Bv/Bv_groups/outgroup_list.txt ~/transcriptomes/baits_Bv/Bv_groups/ none none none
tbaits_intron_removal.py LocusKeyFile BlastFilePre SeqFilePre SeqFolder BlastFolder OutFolder OutFilePre OGFileName AlFolder AlFilePre AlFilePost FilePath
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
tbaits_intron_removal.py version 1.0
This script removes introns from sequences, even when the intron/exon boundaries
may not be in the same places as they are in the reference sequences.
tbaits_intron_removal.py
[tab-delimited file with paralog name (tab) gene family name]
[prefix for blast files, or "none", if none]
[prefix for sequence files, or "none", if none]
[folder where sequence files are found]
[folder where blast output files are found, or "same", if it is the same folder 
as where the sequences are found]
[folder where output files should be written, or "same", if it is the same 
folder as where the sequences are found]
[prefix for output files, or "none", if none]
[for phylogenetic analysis: file with the list of outgroups for each locus]
[folder containing the template alignments and backbone trees]
[alignment file prefix or "none"]
[alignment file ending or "none"]
[path to the fasta_to_phylip.py script or "none" if it is in the default path]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 13:
	sys.exit("ERROR!  This script requires 12 additional arguments, and you supplied %d.  %s\n" % (len(sys.argv)-1, Usage))
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
OGFileName = sys.argv[8]
AlFolder = sys.argv[9]
AlFilePre = sys.argv[10]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[11]
if AlFilePost == "none":
	AlFilePost = ""
FilePath = sys.argv[12]
if FilePath[-1] != "/":
	FilePath += "/"
if FilePath == "none/":
	FilePath = ""

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
	print ("%d exons, each of which had at least %d blast hits, were found in the various contigs.\n" % (NumExons, MinNumHits))
	sys.stderr.write("%d exons, each of which had at least %d blast hits, were found in the various contigs.\n" % (NumExons, MinNumHits))
	return TempDict

#IntronPruner takes a dictionary of sequences grouped by locus (output by LocusSeqGetter),
#the information about whether a sequence is in the same direction as the baits or 
#reversed (output by BlastHitsCounter), and information on exon position (output by
#ExonFinder) and makes new sequences with each intron replaced by two dashes.
#A separate dictionary is made for sequences that are sometimes match in a forward
#direction and sometimes match in a reverse direction.  These two dictionaries are
#combined into one for easier printing.
def IntronPruner(SeqDict, FRDict, PosDict):
	TempDict = defaultdict(dict)
	NonAmbigDict = defaultdict(dict)
	AmbigDict = defaultdict(dict)
	NumSeqs = 0
	NumExons = 0
	AmbigSeqs = 0
	NoExons = 0
	for Locus in SeqDict:
		for Contig in SeqDict[Locus]:
			OldSeq = SeqDict[Locus][Contig]
			NewSeq = ""
			#print ("%s: %s" % (Locus,Contig))
			for ExonTuple in PosDict[Locus][Contig]:
				EStart = ExonTuple[0]-1
				EEnd = ExonTuple[1]
				NewSeq += OldSeq[EStart:EEnd]+"--"
				NumExons += 1
			if NewSeq != "":
				NewSeq = NewSeq[:-2]
				if FRDict[Locus][Contig]['F'] == 0: #the contig is definitely reversed
					NewSeq2 = Seq(NewSeq, IUPAC.unambiguous_dna)
					NewSeq = str(NewSeq2.reverse_complement())
					NonAmbigDict[Locus][Contig] = NewSeq
					NumSeqs += 1
				elif FRDict[Locus][Contig]['R'] == 0: # the contig is definitely not reversed
					NonAmbigDict[Locus][Contig] = NewSeq
					NumSeqs += 1
				else: #it is not clear
					print("WARNING!! There is some ambiguity for contig %s from locus %s, because it matched in the reverse direction %d times and in the forward direction %d times.\n" % (Contig, Locus, FRDict[Locus][Contig]['R'], FRDict[Locus][Contig]['F']))
					sys.stderr.write("WARNING!! There is some ambiguity for contig %s from locus %s, because it matched in the reverse direction %d times and in the forward direction %d times.\n" % (Contig, Locus, FRDict[Locus][Contig]['R'], FRDict[Locus][Contig]['F']))
					NewName = "F"+str(FRDict[Locus][Contig]['F'])+"_R"+str(FRDict[Locus][Contig]['R'])+"_"+Contig
					AmbigDict[Locus][NewName] = NewSeq
					AmbigSeqs += 1				
			else:
				NoExons += 1
	TempDict[''] = NonAmbigDict
	TempDict['Ambig'] = AmbigDict
	print("%d exons in %d unambiguous contigs and %d ambiguous contigs were transcribed.\n" % (NumExons, NumSeqs, AmbigSeqs))
	sys.stderr.write("%d exons in %d unambiguous contigs and %d ambiguous contigs were transcribed.\n" % (NumExons, NumSeqs, AmbigSeqs))
	print("%d contigs did not blast strongly to the exon database.\n" % (NoExons))
	sys.stderr.write("%d contigs did not blast strongly to the exon database.\n" % (NoExons))
	return TempDict

#LocusGroupSeqWriter writes sequences in the desired format to files from a
#dictionary where the top-level key is the group, the second-level key is the locus,
#the third-level key is the sequence name, and the value is the sequence.
#The top-level key can simply be "".
def LocusGroupSeqWriter(SeqDict, Folder, Prefix, Suffix, SeqSuffix, SeqFormat):
	TempDict = defaultdict(dict)
	for LocusGroup in SeqDict:
		NumFiles = 0
		for Locus in SeqDict[LocusGroup]:
			OutFileName = Folder+Prefix+LocusGroup+Locus+Suffix
			TempDict[LocusGroup][Locus] = Prefix+LocusGroup+Locus
			OutFile = open(OutFileName, 'w')
			for Contig in SeqDict[LocusGroup][Locus]:
				Record1 = SeqRecord(seq=Seq(SeqDict[LocusGroup][Locus][Contig], IUPAC), id = Contig+SeqSuffix, description = "")
				SeqIO.write(Record1, OutFile, SeqFormat)
			OutFile.close()
			NumFiles += 1
		print("%d sequence files were written for the group %s, with names such as %s.\n" % (NumFiles, LocusGroup, OutFileName))
		sys.stderr.write("%d sequence files were written for the group %s, with names such as %s.\n" % (NumFiles, LocusGroup, OutFileName))
	return TempDict

#SeqStatsWriter reads the sequences and calculates some statistics.
#original?
def SeqStatsWriter(SeqDict, Folder, Pre, Delimitter):
	TempDict = defaultdict(dict) #TempDict[Ind][Locus] = NumContigs
	MeanTempDict = { } #MeanTempDict[Locus][Ind] = [list of contig lengths]
	#making the dictionary showing the number of contigs per locus
	for Locus in SeqDict['']:
		MeanTempDict[Locus] = defaultdict(list)
		for ContigName in SeqDict[''][Locus]:
			IndName = ContigName.split(Delimitter)[0]
			MeanTempDict[Locus][IndName].append(SeqLenFinder(SeqDict[''][Locus][ContigName]))
			try:
				TempDict[IndName][Locus] += 1
			except KeyError:
				TempDict[IndName][Locus] = 1
	#writing the dictionary to the file (TempDict)
	TempIndList = sorted(TempDict.keys())
	OutList = ['Locus\t'+"\t".join(TempIndList)]
	for Locus in SeqDict['']:
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
	for Locus in SeqDict['']:
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
		
#MRScriptWriter writes an output script for analysis of the sequences using mafft and raxml
def MRScriptWriter(SeqFileDict, Folder, Prefix, AFolder, APre, APost, Path):
	for LocusGroup in SeqFileDict:
		OutList = [ ]
		OutLocusList = [ ]
		if LocusGroup != "Ambig":
			for Locus in SeqFileDict[LocusGroup]:
				NamePart = SeqFileDict[LocusGroup][Locus]
				#Line = "rm "+Folder+NamePart+"_exons_al.fa && "
				#Line += "rm "+Folder+"RAxML*"+NamePart+"\n"
				#OutList.append(Line)
				Line = "mafft --addfragments "+Folder+NamePart+".fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+NamePart+"_exons_al.fa && "
				Line += Path+"fasta_to_phylip.py "+Folder+NamePart+"_exons_al.fa && "
				Line += "raxmlHPC -f v -s "+Folder+NamePart+"_exons_al.phy -n "+NamePart+" -t "+AFolder+"RAxML_bipartitions."+APre+Locus+" -m GTRCAT -o "+OGDict[Locus]+" -w "+Folder+"\n"
				OutList.append(Line)
				OutLocusList.append(Locus+"\n")
			OutFileName = Folder+Prefix+LocusGroup+"Analysis_Script.sh"
			OutFileWriting(OutFileName, OutList)
			print("The shell script for analyzing this %s group of sequences data was written to %s.\n" % (LocusGroup, OutFileName))
			sys.stderr.write("The shell script for analyzing this %s group of sequences data was written to %s.\n" % (LocusGroup, OutFileName))
			OutFileName = Folder+Prefix+LocusGroup+"Locus_List.txt"
			OutFileWriting(OutFileName,OutLocusList)
			print("The list of loci was written to the file %s.\n" % (OutFileName))
			sys.stderr.write("The list of loci was written to the file %s.\n" % (OutFileName))
		else:
			for Locus in SeqFileDict[LocusGroup]:
				Line = Locus+"\t"+SeqFileDict[LocusGroup][Locus]+"\n"
				OutList.append(Line)
			OutFileName = Folder+Prefix+LocusGroup+"_Ambiguous_Contigs.txt"
			OutFileWriting(OutFileName, OutList)
			print("The list of files of ambiguous sequences was written to %s.\n" % (OutFileName))
			sys.stderr.write("The list of files of ambiguous sequences was written to %s.\n" % (OutFileName))


#Read the BlastFileList and the LocusKeyFile and make dictionaries from them:
LocusDict = DictFromFile(LocusKeyFile)
LocusList = ListFromDictValues(LocusDict)
OGDict = DictFromFile(OGFileName)
SeqFileList = [SeqFilePre+Locus+".fa" for Locus in LocusList]
BlastFileList = [BlastFilePre+Locus+".out" for Locus in LocusList]

#read contig files, saving the contigs, and preparing the dictionary that will tell how many times blast found a given part of a sequence

ContigDict = LocusSeqGetter(SeqFileList, SeqFolder, SeqFilePre, ".fa", 'fasta')
SeqPosDicta = SeqRangeDictMaker(ContigDict)

#read blast files and figure out how many times each base in the contig file was hit

SeqPosDict = BlastHitsCounter(SeqPosDicta, BlastFileList, BlastFolder, 1e-12, LocusDict)

#determine which positions had enough hits to be exons

ExonPosDict = ExonFinder(SeqPosDict, 2)

#remove the exons and decide what to do with them when we reach that point

OutDict = IntronPruner(ContigDict, SeqPosDict, ExonPosDict)

#Figure out which loci had sequences
SeqStatsWriter(OutDict, OutFolder, OutFilePre, "-")

#Write sequence files and shell scripts

FileDict = LocusGroupSeqWriter(OutDict, OutFolder, OutFilePre, '.fa', "_exons", 'fasta')
MRScriptWriter(FileDict, OutFolder, OutFilePre, AlFolder, AlFilePre, AlFilePost, FilePath)
