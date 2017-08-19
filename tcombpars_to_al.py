#! /usr/bin/env python

#tcombpars_to_al.py version 1.0 4 Aug. 2017 Abby Moore
#This script takes the alignments produced by tparcomb_final.py or tambig_seq_renaming.py and
#makes alignments of the subset of the sequences that are present in either the species tree
#or a list of taxa.  It also removes portions that are mainly gaps. 
#Then it writes a shell script to align the sequence files.
#Modified from tcombpars_to_trees.py

import sys
from collections import defaultdict
import dendropy
from Bio.Seq import Seq #to edit sequences
from Bio import AlignIO #to read and parse sequence alignments
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment #to make multiple sequence alignments
import subprocess #We want to be able to talk to the command line.

'''
tcombpars_to_trees.py LocusListFN SpeciesTreeFN AlFolder1 AlFilePre1 AlFilePost1 AlFolder2 AlFilePre2 AlFilePost2 OutFolder OutFilePre OutFilePost ScriptPath Mode[Parallel, Array] NextScriptFN(if Mode == Array)   
'''

Usage = '''
tcombpars_to_al.py version 1.0
This script takes the set of sequences we want to analyze further from gene
family alignments.  It writes a script to align those gene families.
tcombpars_to_al.py
[file containing the list of loci]
[the species tree]
[alignment folder--This one will be tried first and the second one will
not be tried if an alignment is found here!]
[prefix for first set of alignments, or "none", if none]
[suffix for first set of alignments, or "none", if none]
[backbone alignment folder, or "none", if none]
[prefix for backbone alignments, or "none", if none]
[suffix for backbone alignments, or "none", if none]
[output folder, or "same" if the same as the alignment folder]
[prefix for the output files, or "none", if none]
[suffix for the output files, or "none", if none]
[path to the scripts, or "none", if they are on the default path]
[mode: Parallel or Array]
[way to find individuals: Tree or List]
[method for finding gaps: "Empty" if only gaps count, "List" if all positions 
that only include the individuals in the list should be removed, an integer
if all positions with that many or fewer individuals should be removed, or a 
proportion (< 1) if all positions with that proportion of the individuals or
fewer should be removed]
[file name for the outgroup dictionary]
[file name for the subsequent script, only if the mode is Array, or "none", if none]
'''

print("%s\n" % (" ".join(sys.argv)))
ModeList = ['Parallel', 'Array']
IModeList = ['Tree', 'List']

if len(sys.argv) < 17:
	sys.exit("ERROR!  This script requires 16 or 17 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
LocusListFN = sys.argv[1]
SpeciesTreeFN = sys.argv[2]
AlFolder1 = sys.argv[3]
if AlFolder1[-1] != "/":
	AlFolder1 += "/"
AlFilePre1 = sys.argv[4]
if AlFilePre1 == "none":
	AlFilePre1 = ""
AlFilePost1 = sys.argv[5]
if AlFilePost1 == "none":
	AlFilePost1 = ""
BBAlFolder = sys.argv[6]
if BBAlFolder[-1] != "/":
	BBAlFolder += "/"
BBAlFilePre = sys.argv[7]
if BBAlFilePre == "none":
	BBAlFilePre = ""
BBAlFilePost = sys.argv[8]
if BBAlFilePost == "none":
	BBAlFilePost = ""
OutFolder = sys.argv[9]
if OutFolder == "same":
	OutFolder = AlFolder1
elif OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[10]
if OutFilePre == "none":
	OutFilePre = ""
OutFilePost = sys.argv[11]
if OutFilePost == "none":
	OutFilePost = ""
ScriptPath = sys.argv[12]
if ScriptPath[-1] != "/":
	ScriptPath += "/"
if ScriptPath == "none/":
	ScriptPath = ""
Mode = sys.argv[13]
if Mode not in ModeList:
	sys.exit("ERROR!!  You wanted the mode %s, but the mode must be one of the following: %s.\n%s" % (Mode, ", ".join(ModeList), Usage))
IMode = sys.argv[14]
if IMode not in IModeList:
	sys.exit("ERROR!!  You wanted the mode %s, but the mode must be one of the following: %s.\n%s" % (IMode, ", ".join(IModeList), Usage))
Method = sys.argv[15]
OGDictFN = sys.argv[16]
if Mode == "Array":
	NextScriptFN = sys.argv[17]

Vociferous = False

######################################################################################

#SeqFileReading reads a sequence file and puts the sequences in a dictionary.
#not original, but not sure where this is from
def SeqFileReading(FileName, SeqFormat):
	DictTemp = { }#DictTemp[SeqName] = Seq
	InFile = open(FileName, 'rU')
	for record in SeqIO.parse(InFile, SeqFormat):
		DictTemp[record.id] = str(record.seq)
	InFile.close()
	return DictTemp
	#This is various sequence dictionaries

#SeqFileWriting writes sequence files from dictionaries.
#not original, but not sure where this is from
def SeqFileWriting(FileName, SDict, SeqFormat):
	OutFile = open(FileName, 'w')
	for SeqName in SDict:
		Record1 = SeqRecord(seq=Seq(SDict[SeqName], IUPAC), id = SeqName, description = "")
		SeqIO.write(Record1, OutFile, SeqFormat)
	OutFile.close()
	return

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
#a modified version of the function in from tbaits_intron_removal.py from tcontig_selection.py
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
	#This is OutGroupDict

#IndListfromTree makes a list of the sequences present in a tree
#This is no longer being used, but I will leave it for now, since this is the only script it is in.
#original
def IndListfromTree(FileName):
	InFile = open(FileName, 'rU')
	for Line in InFile:
		TempTree = dendropy.Tree.get_from_string(Line.strip('\n').strip('\r'), schema = 'newick', preserve_underscores=True)
	InFile.close()
	ListTemp = [Node.taxon.label for Node in TempTree.leaf_nodes()]
	return ListTemp
	#This is IndList.

#GetIndName gets the name of the individual from a sequence, which can either be a contig or a full sequence
#tnotung_homolog_parsing.py
def GetIndName(NodeName):
	#for full sequences:
	if (len(str(NodeName).split(".")) > 2):
		IndNameTemp = str(NodeName).split(".")[0]
	#for contigs
	else:
		IndNameTemp = str(NodeName).split("-")[0]
	return IndNameTemp
	#This is IndName

#OGDictfromTree makes an ordered list of individuals in a tree, which can be used to choose outgroup sequences
#original
def OGDictfromTree(FileName):
	DictTemp = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		TempTree = dendropy.Tree.get_from_string(Line.strip('\n').strip('\r'), schema = 'newick', preserve_underscores=True)
	InFile.close()
	#print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	IndNum = 0
	#if we are running on Dendropy 4.xx
	try:
		for Node in TempTree.levelorder_node_iter():
			if Node.is_leaf():
				IndNum += 1
				DictTemp[Node.taxon.label] = IndNum
	#if we are running on Dendropy 3.xx
	except AttributeError:
		for Node in TempTree.level_order_node_iter():
			if Node.is_leaf():
				IndNum += 1
				DictTemp[Node.taxon.label] = IndNum
	
	return (DictTemp, IndNum)
	#DictTemp is OGIndDict and IndNum is NumInds

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

#GapFinder finds the gaps in an alignment.  It returns a dictionary in which the gaps are numbered consecutively
#(starting with 0) of the following form: DictTemp[GapNum]['start'/'end']: position in the alignment where the gap starts
#or ends, respectively
#Meth is the method, "empty" for only empty positions removed, "list" for remove positions in which members of the list are the only
#ones present, and a number for remove positions in which that many or fewer individuals are present
#ListTemp is the list of individuals to remove, or a blank list.
#from tal_combiner.py
def GapFinder(Meth, AlignTemp, ListTemp):
	#look through each position of the alignment
	DictTemp = defaultdict(dict)#DictTemp[GapNum]['start'/'end'] = SeqPos
	AlLength = len(AlignTemp)
	if Meth not in ["Empty", "List"]:
		#if you want at least a certain number of individuals to be present
		if float(Meth) > 1:
			MinIndsPresent = int(Meth)
		#if you want at least a certain proportion of individuals to be present
		else:
			MinIndsPresent = float(Meth)*AlLength
	GapNum = 0
	GapStart = 0
	InGap = True
	GapLength = 0
	for SeqPos in range(0, AlignTemp.get_alignment_length()):
		#determine whether or not there are nucleotides at that position
		GapOnly = True
		if Meth == "Empty":
			for record in AlignTemp:
				if (record[SeqPos] != '-') and (record[SeqPos] != '?'):
					GapOnly = False
		elif Meth == "List":
			for record in AlignTemp:
				if record.id not in ListTemp:
					if (record[SeqPos] != '-') and (record[SeqPos] != '?'):
						GapOnly = False
		else:
			NumInds = 0
			for record in AlignTemp:
				if (record[SeqPos] != '-') and (record[SeqPos] != '?'):
					NumInds += 1
			if NumInds > MinIndsPresent:
				GapOnly = False 
		#four possible combinations of currently in a gap or not and current site is a gap or not
		#if this is a start of a new sequence segment
		if (GapOnly == False) and (InGap == True):
			DictTemp[GapNum]['start'] = GapStart
			DictTemp[GapNum]['end'] = SeqPos-1
			GapLength += (SeqPos-GapStart)
			GapNum += 1
			InGap = False
		#the start of a new gap
		elif (GapOnly == True) and (InGap == False):
			GapStart = SeqPos
			InGap = True
		#We do not need to deal with the other two options:
		#inA the middle of a segment [(GapOnly == False) and (InGap == False)]
		#or in the middle of a gap [(GapOnly == True) and (InGap == True)]
	#check to see if there is a gap at the end of the alignment
	if (InGap == True):
		DictTemp[GapNum]['start'] = GapStart
		DictTemp[GapNum]['end'] = SeqPos
		GapLength += (SeqPos-GapStart+1)
	if Vociferous == True: print("%d gaps were found, with a total length of %d.\n" % (GapNum+1, GapLength))
	return DictTemp
	#This is GapDict

#AlignmentGapRemoving takes an alignment and removes all of the characters that are gaps
#in all sequences.  It uses the results of GapFinder.
#from tal_combiner.py
def AlignmentGapRemoving(AlignTemp, DictTemp):
	AlEnd = AlignTemp.get_alignment_length()-1
	if Vociferous == True: print("Before removing any gaps, the alignment is %d bases long.\n" % (AlEnd+1))
	for GapNum in sorted(DictTemp.keys(), reverse=True):
		#print("Before removing gap %d: gap length: %d, alignment length %d\n" % (GapNum, DictTemp[GapNum]['end']-DictTemp[GapNum]['start']+1, AlignTemp.get_alignment_length()))
		#if this gap encompasses the end of the alignment, cut that bit off of the end of the alignment
		if DictTemp[GapNum]['end'] == AlEnd:
			AlignTemp = AlignTemp[:, :DictTemp[GapNum]['start']]
		#if this gap encompasses the start of the alignment, cut off the start
		elif DictTemp[GapNum]['start'] == 0:
			AlignTemp = AlignTemp[:, (DictTemp[GapNum]['end']+1):]
		#for all gaps in the middle of the alignment
		else:
			AlignTemp = AlignTemp[:, :DictTemp[GapNum]['start']] + AlignTemp[:, (DictTemp[GapNum]['end']+1):]
		#print("After removing gap %d: alignment length %d\n" % (GapNum, AlignTemp.get_alignment_length()))
	if Vociferous == True: print("After removing %d gaps, the alignment is now %d bases long.\n" % (len(DictTemp), AlignTemp.get_alignment_length()))
	return AlignTemp
	#This is now PAlign

#MissingIndRemoving goes through an alignment once the gaps have been removed and takes out any
#individuals that then consist entirely of missing data
#from tal_combiner.py
def MissingIndRemoving(AlignTemp):
	MissingSeqDict = { }
	ListTemp = [ ]
	#Go through the entire alignment
	SeqNum = 0
	for Record in AlignTemp:
		#If any of the rows consist entirely of missing data,
		if list(set(list(str(Record.seq)))) == ['-']:
			#add that sequence to the list of missing sequences
			MissingSeqDict[SeqNum] = Record.id
		#if not, add it to the list of good individuals
		else:
			ListTemp.append(Record.id)
		SeqNum += 1
	#If there are missing sequences, remove them from the alignment
	if MissingSeqDict != { }:
		MissingSeqList = sorted(MissingSeqDict.keys())
		if len(MissingSeqList) == 1:
			NewAlign = AlignTemp[0:MissingSeqList[0]]
			NewAlign.extend(AlignTemp[(MissingSeqList[0]+1):(len(AlignTemp))])
			AlignTemp = NewAlign
		elif len(MissingSeqList) > 1:
			NewAlign = AlignTemp[0:MissingSeqList[0]]
			for PosNum in list(range(len(MissingSeqList))):
				if PosNum != len(MissingSeqList)-1:
					NewAlign.extend(AlignTemp[(MissingSeqList[PosNum]+1):MissingSeqList[PosNum+1]])
				else:
					NewAlign.extend(AlignTemp[(MissingSeqList[PosNum]+1):len(AlignTemp)])
			AlignTemp = NewAlign
	return (AlignTemp, ListTemp)			

##############################################################################################################################

#reading the list of loci
LocusList = CaptureColumn(LocusListFN, 0)
print("The list of %d loci was read from file %s.\n" % (len(LocusList), LocusListFN))
sys.stderr.write("The list of %d loci was read from file %s.\n" % (len(LocusList), LocusListFN))

OutGroupDict = DictFromFile(OGDictFN, 0, 1)


#reading the list of individuals from the tree or from a list
if IMode == "Tree":
	IndList = IndListfromTree(SpeciesTreeFN)
	print("The names of %d individuals were read from the tree %s.\n" % (len(IndList), SpeciesTreeFN))
	sys.stderr.write("The names of %d individuals were read from the tree %s.\n" % (len(IndList), SpeciesTreeFN))
elif IMode == "List":
	IndList = CaptureColumn(SpeciesTreeFN,0)

#currently we don't have the option of using a list of long sequences to find gaps
LongSeqList = [ ]
	

#reading the old sequence files and writing new ones that just have the individuals we want to analyze further
LocusListRemove = [ ]
NumSeqsDict = defaultdict(list)
SeqsperLocus = { }
for Locus in LocusList:
	FileExists = True
	InFileName1 = AlFolder1+AlFilePre1+Locus+AlFilePost1
	try:
		SeqDictIn = SeqFileReading(InFileName1, 'fasta')
		if len(SeqDictIn) == 0:
			FileExists = False
		if Vociferous == True: print("%d sequences were read for locus %s from file %s.\n" % (len(SeqDictIn), Locus, InFileName1))
	except IOError:
		LocusListRemove.append(Locus)
		print("Locus %s removed from the list." % (Locus))
		FileExists = False
	if FileExists:
		SeqDictOut = { }
		for SeqName in SeqDictIn:
			if GetIndName(SeqName) in IndList:
				SeqDictOut[SeqName] = SeqDictIn[SeqName]
		if Vociferous == True:
			print("%d of these sequences belong to the individuals of interest.\n" % (len(SeqDictOut)))
		NumSeqsDict[len(SeqDictOut)].append(Locus)
		SeqsperLocus[Locus] = len(SeqDictOut)
		OutFileName = OutFolder+OutFilePre+Locus+OutFilePost+".fa"
		AlignGaps = 0
		for SeqName in SeqDictOut:
			Record1 = SeqRecord(seq=Seq(SeqDictOut[SeqName]), id = SeqName, description = "")
			try:
				AlignGaps.append(Record1)
			except AttributeError:
				AlignGaps = MultipleSeqAlignment([Record1])
		GapPosDict = GapFinder(Method, AlignGaps, LongSeqList)
		AlignNoGaps = AlignmentGapRemoving(AlignGaps, GapPosDict)
		(AlignOut, SeqsUsing) = MissingIndRemoving(AlignNoGaps)
		if len(AlignNoGaps) != len(AlignOut):
			print("%d sequences that were only gaps were removed." % (len(AlignNoGaps)-len(AlignOut)))
		AlignIO.write(AlignOut, OutFileName, "fasta")
print("New sequence files for %d loci were written, with names such as %s.\n" % (len(LocusList), OutFileName))
sys.stderr.write("New sequence files for %d loci were written, with names such as %s.\n" % (len(LocusList), OutFileName))

#re-ordering the LocusList so the files are in order from largest to smallest
LocusList = [ ]
for NumSeqs in sorted(NumSeqsDict.keys(), reverse=True):
	if NumSeqs != 0:
		LocusList += NumSeqsDict[NumSeqs]

#writing the script to analyze the data further
if Mode == 'Parallel':
	Script1FileName = OutFolder+OutFilePre+"analysis_script1.sh"
	Script2FileName = OutFolder+OutFilePre+"analysis_script2.sh"
	OutList1 = ["#! /bin/bash\n\n"]
	OutList2 = []
	for Locus in LocusList:
		#Add second alignment
		if BBAlFolder != "none/":
			Line = ScriptPath+"alignment_combining.py "+BBAlFolder+BBAlFilePre+Locus+BBAlFilePost+".fa fasta "+OutFolder+OutFilePre+Locus+OutFilePost+".fa fasta "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa fasta\n"
			OutList2.append(Line)
		#Align
		Line = "rm "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa "+OutFolder+"RAxML_*."+OutFilePre+Locus+OutFilePost+"\n"
		OutList1.append(Line)
		Line = "mafft --localpair --maxiterate 1000 --quiet "+OutFolder+OutFilePre+Locus+OutFilePost+".fa > "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa && raxmlHPC -f a -s "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.phy -n "+OutFilePre+Locus+OutFilePost+" -m GTRCAT -p 1234 -N 100 -x 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
		OutList2.append(Line)
	Line = "cat "+Script2FileName+" | parallel --joblog "+OutFolder+OutFilePre+"parallel_log2.log\n"
	OutList1.append(Line)
	OutFileWriting(Script1FileName, OutList1)
	OutFileWriting(Script2FileName, OutList2)
elif Mode == 'Array':
	SBatchList = [ ]
	for Locus in LocusList:
		OutList = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+Locus+"\n"]
		if SeqsperLocus[Locus] < 100:
			Line = "#SBATCH -t 2:00:00\n#--mem=4GB\n"
		elif SeqsperLocus[Locus] < 200:
			Line = "#SBATCH -t 3:00:00\n#--mem=8GB\n"
		else:
			Line = "#SBATCH -t 4:00:00\n#--mem=16GB\n"
		Line += "module load mafft\nmodule load raxml\n"
		OutList.append(Line)
		#Add second alignment
		if BBAlFolder != "none/":
			Line = ScriptPath+"alignment_combining.py "+BBAlFolder+BBAlFilePre+Locus+BBAlFilePost+".fa fasta "+OutFolder+OutFilePre+Locus+OutFilePost+".fa fasta "+OutFolder+OutFilePre+Locus+OutFilePost+".fa fasta\n"
			OutList.append(Line)
		#Align
		Line = "rm "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa "+OutFolder+"RAxML_*."+OutFilePre+Locus+OutFilePost+"\n"
		OutList.append(Line)
		Line = "mafft --localpair --maxiterate 1000 --quiet "+OutFolder+OutFilePre+Locus+OutFilePost+".fa > "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa\n"
		OutList.append(Line)
		Line ="raxmlHPC -f v -s "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa -n "+OutFilePre+Locus+OutFilePost+" -m GTRCAT -t "+BBAlFolder+"RAxML_bipartitions."+BBAlFilePre+Locus+" -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
		OutList.append(Line)
		OutFileName = OutFolder+OutFilePre+Locus+"script.sh"
		OutFileWriting(OutFileName, OutList)
		OutLine = "sbatch "+OutFileName
		SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		SBatchList.append(SBatchOut.strip('\r').strip('\n').split(" ")[-1])
	print("%d separate scripts for alignment, one for each locus, were submitted.\n" % (len(SBatchList)))
	sys.stderr.write("%d separate scripts for alignment, one for each locus, were submitted.\n" % (len(SBatchList)))
	if NextScriptFN != "none":
		OutLine = "sbatch -d afterok:"+":".join(SBatchList)+" "+NextScriptFN
		SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
		print("The subsequent script, %s, will be run with the job id %s, when the previous scripts have finished.\n" % (NextScriptFN, JobIDPrev))
		sys.stderr.write("The subsequent script, %s, will be run with the job id %s, when the previous scripts have finished.\n" % (NextScriptFN, JobIDPrev))
