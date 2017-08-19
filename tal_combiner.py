#! /usr/bin/env python

#tal_combiner.py version 2.0 2 Nov 2015 Abby Moore
#version 2.0 removes the gaps from the alignments before writing them to files
#It reads the locus_list_out.txt file from tnotung_homolog_parsing.py

#tal_combiner.py ~/transcriptomes/Anacampserotaceae/Anac_seqs_final/LocusList.txt ~/transcriptomes/Anacampserotaceae/Anac_seqs_final/ Anac_ _OGa.fa fasta fasta Pereskia ~/transcriptomes/Anacampserotaceae/Anac_seqs_final/seqs_pruned/ pAnac_ separate 0 0 all
#tal_combiner.py LocusListFileName ScriptFolder SeqFolder SeqFilePre SeqFilePost SeqFormatIn SeqFormatOut OGName OutFolder OutFilePre OutMode[combined, separate, separatemb, combinednobs] MaxIndsMissing MaxSeqsMissing IndsUsingFileName[all] Method[Empty, List, integer] Mode[Array, Parallel] OGTreeFN[if OGName == tree] LongSeqListFN[if Method == List]

Usage = '''
tal_combiner renames the sequences in the alignments so that they are named
after the individual only instead of having all of the other information.
It identifies the individual name by taking everything before the first period.
It can either output several separate alignments with the same names or it can
output one concatenated alignment.
tal_combiner.py
[file name for the list of loci to be combined]
[folder where the scripts are found, or none, if they are in the path]
[folder where the separate alignments for each locus are found]
[prefix for alignment files, or none, if none]
[suffix for alignment files, or none, if none]
[format of alignment files--anything biopython can read]
[desired format for output alignments--anything biopython can write]
[name of the outgroup; or none, if none; or tree, if the outgroup should be the
basal-most individual that is present in the alignment, based on a species tree
--outgroup must be present in all alignments, if it is named specifically]
[folder to which the new alingments should be output, or same, if same as the
input folder]
[desired prefix for new alignments, or none, if none]
[mode of output: combined, separate, combinednobs, if I am just making an
intermediate species tree, or separatemb, if I want to run MrBayes for bucky,
in which case we need to have the missing taxa included as question marks]
[the maximum number of missing individuals a locus can have and still be used]
[the maximum number of missing loci an individual can have and still be used]
[file name for a list of which individuals should be used, or all, if all 
individuals should be used]
[method for finding gaps: "Empty" if only gaps count, "List" if all positions 
that only include the individuals in the list should be removed, an integer
if all positions with that many or fewer individuals should be removed, or a 
proportion (< 1) if all positions with that proportion of the individuals or
fewer should be removed]
[mode, either Array or Parallel]
[tree file for outgroup, if OGName is tree]
[list of individuals whose sequences should be removed if they are the only ones
present at a given site, if Method is List]
'''

import sys
import subprocess
from collections import defaultdict
from Bio import AlignIO #to read and parse sequence alignments
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment
import dendropy #to get outgroups from trees
import copy#to prune trees without destroying the original

OutModeList = ['combined', 'separate', 'separatemb', 'combinednobs']
ModeList = ['Array', 'Parallel']

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 17:
	sys.exit("ERROR!  This script requires at least 16 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
LocusListFileName = sys.argv[1]
ScriptFolder = sys.argv[2]
if ScriptFolder[-1] != "/":
	ScriptFolder += "/"
if ScriptFolder == "none/":
	ScriptFolder = ""
SeqFolder = sys.argv[3]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
SeqFilePre = sys.argv[4]
if SeqFilePre == "none":
	SeqFilePre = ""
SeqFilePost = sys.argv[5]
if SeqFilePost == "none":
	SeqFilePost = ""
SeqFormatIn = sys.argv[6]
SeqFormatOut = sys.argv[7]
OGName = sys.argv[8]
OutFolder = sys.argv[9]
if OutFolder == "same":
	OutFolder = SeqFolder
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[10]
if OutFilePre == "none":
	OutFilePre = ""
OutMode = sys.argv[11]
if (OutMode in OutModeList) == False:
	sys.exit("ERROR!  For output mode, you wrote %s, but it can only be %s.\n%s" % (OutMode, " ".join(OutModeList), Usage))
MaxIndsMissing = int(sys.argv[12])
MaxSeqsMissing = int(sys.argv[13])
IndsUsingFileName = sys.argv[14]
Method = sys.argv[15]
Mode = sys.argv[16]
if Mode not in ModeList:
	sys.exit("ERROR!!  You wanted the mode to be %s, but it must be one of the following: %s.\n%s" % (Mode, ", ".join(ModeList), Usage))
if OGName == "tree":
	OGTreeFN = sys.argv[17]
	if Method == "List":
		LongSeqListFN = sys.argv[18]
elif Method == "List":
	LongSeqListFN = sys.argv[17]


Verbose = False
Verbose = True

###################################################################################

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
	#This is LocusList and IndsUsing

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
	if Verbose == True: print("%d gaps were found, with a total length of %d.\n" % (GapNum+1, GapLength))
	return DictTemp
	#This is GapDict

#AlignmentGapRemoving takes an alignment and removes all of the characters that are gaps
#in all sequences.  It uses the results of GapFinder.
def AlignmentGapRemoving(AlignTemp, DictTemp):
	AlEnd = AlignTemp.get_alignment_length()-1
	if Verbose == True: print("Before removing any gaps, the alignment is %d bases long.\n" % (AlEnd+1))
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
	if Verbose == True: print("After removing %d gaps, the alignment is now %d bases long.\n" % (len(DictTemp), AlignTemp.get_alignment_length()))
	return AlignTemp
	#This is now PAlign

#MissingIndRemoving goes through an alignment once the gaps have been removed and takes out any
#individuals that then consist entirely of missing data
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
		

#LocusAlSeqGetter reads a series of alignments and makes a dictionary of the sequences
#that have been classified according to locus.
def LocusAlSeqGetter(LList,Folder,FilePre,FilePost,SeqFormat):
	TempDict = defaultdict(dict)#TempDict[Locus][ContigName] = ContigSeq
	for Locus in LList:
		FileName = Folder + FilePre + Locus + FilePost
		try:
			MyAlignment = AlignIO.read(FileName, SeqFormat)
			for seq_record in MyAlignment:
				TempDict[Locus][seq_record.id] = str(seq_record.seq)
				#print("%d sequences for locus %s were read from the file %s.\n" % (len(TempDict[Locus].keys()), Locus, FileName))
		#If this alignment doesn't exist, that is because the paralog only had one sequence, so mafft didn't align it.
		except ValueError:
			FileName = ".".join(FileName.split("_al."))
			try:
				MyAlignment = AlignIO.read(FileName, SeqFormat)
				for seq_record in MyAlignment:
					TempDict[Locus][seq_record.id] = str(seq_record.seq)
					#print("%d sequences for locus %s were read from the file %s.\n" % (len(TempDict[Locus].keys()), Locus, FileName))
				if Verbose == True: print("The new FileName is %s, and it has %d sequences." % (FileName, len(MyAlignment)))
			except ValueError: #Then the file actually doesn't exist.
				"do nothing"
	print("%d sequence files were read, with names such as %s.\n" % (len(LList), FileName))
	sys.stderr.write("%d sequence files were read, with names such as %s.\n" % (len(LList), FileName))
	return TempDict
	#This is SeqDict

#SeqLenFinder finds the length of a sequence in an alignment, while ignoring all gaps
#from tnotung_homolog_parsing.py
def SeqLenFinder(Seq):
	Len = 0
	for Base in Seq:
		if Base != "-":
			Len += 1
	return Len

#OGNamefromTree takes a list of individuals in an alignment and returns the list of outgroups to make a tree
#that corresponds as closely as possible to the species tree
#original
def OGNamefromTree(TempTree,TempList):
	#go through the species tree and prune the taxa that aren't in the list
	LeavesPruning = [ ]
	#the tree needs to be copied or it will change the original and _mess_everything_up_!!!!
	TempTree2 = copy.deepcopy(TempTree)
	#print(TempTree2.as_ascii_plot(show_internal_node_labels=True))
	for Node in TempTree2.leaf_nodes():
		if Node.taxon.label not in TempList:
			#print("Pruning taxon %s.\n" % (Node.taxon.label))
			LeavesPruning.append(Node.taxon)
		#else:
		#	print("Not pruning taxon %s.\n" % (Node.taxon.label))
	TempTree2.prune_taxa(LeavesPruning, update_bipartitions=True)
	#print("%d of the taxa will be pruned.\n" % (len(LeavesPruning)))
	#print(TempTree2.as_ascii_plot(show_internal_node_labels=True))
	#find the basal split of the species tree and take the smaller half
	LeafNodeList = [ ]
	RootNode = TempTree2.seed_node
	for Daughter_Node in RootNode.leaf_nodes():
		LeafNodeList.append(Daughter_Node.taxon.label)
	RootDict = defaultdict(list)
	NodeNum = 1
	for Child_Node in RootNode.child_node_iter():
		for Daughter_Node in Child_Node.leaf_nodes():
			RootDict[NodeNum].append(Daughter_Node.taxon.label)
		NodeNum += 1
	SizeSmallerClade = len(TempList)
	ListOG = [ ]
	RootedInds = 0
	for NodeNum in RootDict:
		RootedInds += len(RootDict[NodeNum])
		if len(RootDict[NodeNum]) < SizeSmallerClade:
			SizeSmallerClade = len(RootDict[NodeNum])
			ListOG = RootDict[NodeNum]
	NameTemp = ",".join(ListOG)
	return(ListOG, NameTemp, RootedInds)

##################################################################################

#read the list of loci to combined
LocusList = CaptureColumn(LocusListFileName, 0)

#read the list of individuals to be used
if IndsUsingFileName != "all":
	IndsUsing = CaptureColumn(IndsUsingFileName, 0)

#read the sequence files for those loci
SeqDict = LocusAlSeqGetter(LocusList, SeqFolder, SeqFilePre, SeqFilePost, SeqFormatIn)

#reading the species tree, if we want to use it to find the outgroup(s)
if OGName == "tree":
	InFile = open(OGTreeFN, 'rU')
	for Line in InFile:
		SppTree = dendropy.Tree.get_from_string(Line.strip('\n').strip('\r'), schema = 'newick', preserve_underscores=True)
	InFile.close()

#reading the list of individuals to exclude from gap calculations, if we want to use that method
if Method == "List":
	LongSeqList = CaptureColumn(LongSeqListFN, 0)
else:
	LongSeqList = [ ]
	

#figure out which individuals have sequences for which loci
SeqIndDict = defaultdict(dict)
IndList = [ ]
#and make a dictionary of sequence lengths for each locus
SeqLenDict = { }
#look through the sequences rearranged by locus, then sequence name
for Locus in SeqDict:
	for SeqName in SeqDict[Locus]:
		#identify which individual each sequence belongs to
		IndName = SeqName.split(".")[0]
		if len(SeqName.split(".")) < 3:
			IndName = SeqName.split("-")[0]
			if len(IndName) > len(SeqName)-2:
				IndName = SeqName.split(".")[0]
		#if we want to use all the individuals
		if IndsUsingFileName == "all":
			#make a new dictionary in which the sequences are arranged by locus, then individual name
			SeqIndDict[Locus][IndName] = SeqDict[Locus][SeqName]
			SeqLenDict[Locus] = len(SeqDict[Locus][SeqName])
			#print("%s: %d" % (Locus, SeqLenDict[Locus]))
			#add the individual to the list of individuals
			IndList.append(IndName)
		#if we want to only use a subset of on the individuals:
		else:
			#first check to see if they are on the list
			if IndName in IndsUsing:
				#and then add them to the dictionary
				SeqIndDict[Locus][IndName] = SeqDict[Locus][SeqName]
				SeqLenDict[Locus] = len(SeqDict[Locus][SeqName])
				#print("%s: %d" % (Locus, SeqLenDict[Locus]))
				#add the individual to the list of individuals
				IndList.append(IndName)
#condense the list of individuals by removing duplicates
IndList = list(set(IndList))

#look through the sequences and make sure each individual has something (either an actual sequence or 
#a set of gaps of the right length
OutSeqDict = defaultdict(dict)
#also determine which sequences are missing for each individual
MissingSeqDict = defaultdict(list)
#and figure out how many bases each individual has
IndSeqLenDict = defaultdict(int)
#go through each individual
for Ind in IndList:
	#try to find the sequence for each locus for that individual
	for Locus in SeqIndDict:
		#add it if it exists
		try:
			OutSeqDict[Locus][Ind] = SeqIndDict[Locus][Ind]
			IndSeqLenDict[Ind] += SeqLenFinder(SeqIndDict[Locus][Ind])
		except KeyError:
			#if not, add gaps to the alignment
			#if OutMode == "separatemb":
			#	OutSeqDict[Locus][Ind] = "?"*SeqLenDict[Locus]
			#else:
			OutSeqDict[Locus][Ind] = "-"*SeqLenDict[Locus]
			#and add that locus to the list of missing sequences for that individual
			MissingSeqDict[Ind].append(Locus)
#determine which sequences will actually be used:
GoodInds = [ ]
BadInds = [ ]
MissingIndDict = { }
GoodSeqs = [ ]
BadSeqs = [ ]
#look through the individuals and eliminate those that have too many missing sequences
for Ind in IndList:
	try:
		NumMissing = len(MissingSeqDict[Ind])
		if NumMissing <= MaxSeqsMissing:
			GoodInds.append(Ind)
		else:
			BadInds.append(Ind)
	except KeyError:
		GoodInds.append(Ind)
if Verbose == True:
	print ("The following %d individuals will be used, because they had %d or fewer missing loci: %s.\n" % (len(GoodInds), MaxSeqsMissing, ", ".join(GoodInds)))
	print ("The following %d individuals will not be used: %s.\n" % (len(BadInds), ", ".join(BadInds)))
print("Of the %d individuals, %d will be used, because they had %d or fewer missing loci, and %d will not be used.\n" % (len(IndList), len(GoodInds), MaxSeqsMissing, len(BadInds)))

#determine how many missing good individuals each locus has
for Locus in LocusList:
	#it needs to be -1 in case none of the individuals have that locus
	MissingIndDict[Locus] = -1
for Ind in GoodInds:
	try:
		for Locus in MissingSeqDict[Ind]:
			MissingIndDict[Locus] += 1
	except KeyError:
		"no missing sequences"
#look through the loci and eliminate those with too many missing (good) individuals:
for Locus in MissingIndDict:
	#version for normal sequences
	if (MissingIndDict[Locus] < MaxIndsMissing) and (MissingIndDict[Locus] != -1):
	#version for species tree making, when all individuals will probably have ITS
	#if (MissingIndDict[Locus] < MaxIndsMissing):
		GoodSeqs.append(Locus)
	else:
		BadSeqs.append(Locus)
if Verbose == True:
	print("The following %d sequences will be used, because they had %d or fewer missing individuals: %s.\n" % (len(GoodSeqs), MaxIndsMissing, ", ".join(GoodSeqs)))
	print("The following %d sequences will not be used: %s.\n" % (len(BadSeqs), ", ".join(BadSeqs)))
print("Of the %d total sequences, %d will be used, because they had %d or fewer missing individuals, and %d will not be used.\n" % (len(LocusList), len(GoodSeqs), MaxIndsMissing, len(BadSeqs)))

#write the file showing how many bases each individual had
OutList = ['Ind_Name\tNum_Bases\n']
for Ind in GoodInds:
	Line = Ind+"\t"+str(IndSeqLenDict[Ind])+"\n"
	OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Len_Ind_Seqs.txt"
OutFileWriting(OutFileName, OutList)

#write the combined alignment
if (OutMode == "combined") or (OutMode == "combinednobs"):
	OutFileName = OutFolder+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs.fa"
	for Ind in GoodInds:
		TempSeq = ""
		for Locus in GoodSeqs:
			TempSeq += OutSeqDict[Locus][Ind] 
		Record1 = SeqRecord(seq=Seq(TempSeq), id = Ind, description = "")
		try:
			AlignGaps.append(Record1)
		except NameError:
			AlignGaps = MultipleSeqAlignment([Record1])
	GapPosDict = GapFinder(Method, AlignGaps, LongSeqList)
	AlignOut = AlignmentGapRemoving(AlignGaps, GapPosDict)
	AlignIO.write(AlignOut, OutFileName, "fasta")
	print("Concatenated sequences for %d individuals were written to the file %s.\n" % (len(AlignOut), OutFileName))
	sys.stderr.write("Concatenated sequences for %d individuals were written to the file %s.\n" % (len(AlignOut), OutFileName))
	if Mode == 'Parallel':
		OutList = [ "#! /bin/bash\n\n"]
	elif Mode == 'Array':
		OutList = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+"concat_tree\n#SBATCH -t 72:00:00\n#SBATCH -n 1\n#SBATCH --mem=32GB\nmodule load raxml\n"] 
	Line = ScriptFolder+"fasta_to_phylip.py "+OutFileName+"\n"
	Line += "rm "+OutFolder+"RAxML_*."+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs\n"
	if OutMode == "combined":
		if OGName == "none":
			Line += "raxmlHPC -s "+OutFolder+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs.phy -n "+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs -m GTRCAT -p 1234 -f a -N 100 -x 1234 -w "+OutFolder+"\n"
		else:
			Line += "raxmlHPC -s "+OutFolder+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs.phy -n "+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+OGName+" -w "+OutFolder+"\n"
	elif OutMode == "combinednobs":
		if OGName == "none":
			Line += "raxmlHPC -s "+OutFolder+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs.phy -n "+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs -m GTRCAT -p 1234 -f d -w "+OutFolder+"\n"
		else:
			Line += "raxmlHPC -s "+OutFolder+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs.phy -n "+OutFilePre+"combined_"+str(len(GoodInds))+"inds_"+str(len(GoodSeqs))+"seqs -m GTRCAT -p 1234 -f d -o "+OGName+" -w "+OutFolder+"\n"
	OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+"analysis_script.sh"
	OutFileWriting(OutFileName, OutList)
	if Mode == "Array":
		OutLine = "sbatch "+OutFileName
		SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		print("The script %s was started with run id %s.\n" % (OutFileName, SBatchOut))
		sys.stderr.write("The script %s was started with run id %s.\n" % (OutFileName, SBatchOut))
#or write the separate alignments
elif OutMode == 'separatemb':
	if Mode == 'Parallel':
		OutList1 = [ "#! /bin/bash\n\ncd "+OutFolder+"\n"]
		OutList2 = [ ]
		for Locus in GoodSeqs:
			#***Probably need a filter for these loci as well***
			OutFileName = OutFolder+OutFilePre+Locus+"_mb.fa"
			AlignGaps = 0
			for Ind in GoodInds:
				#The version if you want to add blank lines for the missing individuals
				'''
				Record1 = SeqRecord(seq=Seq(OutSeqDict[Locus][Ind]), id = Ind, description = "")
				try:
					AlignGaps.append(Record1)
				except AttributeError:
					AlignGaps = MultipleSeqAlignment([Record1])
				'''
				#The version if you want to leave them out:
				try:
					Record1 = SeqRecord(seq=Seq(SeqIndDict[Locus][Ind]), id = Ind, description = "")
					try:
						AlignGaps.append(Record1)
					except AttributeError:
						AlignGaps = MultipleSeqAlignment([Record1])
				except KeyError:
					"do not add that sequence, because it doesn't exist"
			GapPosDict = GapFinder(Method, AlignGaps, LongSeqList)
			AlignNoGaps = AlignmentGapRemoving(AlignGaps, GapPosDict)
			(AlignOut, SeqsUsing) = MissingIndRemoving(AlignNoGaps)
			if len(AlignNoGaps) != len(AlignOut):
				print("%d sequences that were only gaps were removed." % (len(AlignNoGaps)-len(AlignOut)))
			AlignIO.write(AlignOut, OutFileName, "fasta")
			Line = ScriptFolder+"fasta_to_nexus_mb.py "+OutFilePre+Locus+"_mb.fa && "
			Line += "mb < "+OutFilePre+Locus+"_mb_mb.nex > "+OutFilePre+Locus+"_mb.log\n"
			OutList2.append(Line)
		print("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		sys.stderr.write("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		OutFileName1 = OutFolder+OutFilePre+"analysis_script_mb1.sh"
		OutFileName2 = OutFolder+OutFilePre+"analysis_script_mb2.sh"
		Line = "cat "+OutFileName2+" | parallel --joblog "+OutFolder+OutFilePre+"separate_trees.log\n"
		OutList1.append(Line)
		OutFileWriting(OutFileName1, OutList1)
		OutFileWriting(OutFileName2, OutList2)
	elif Mode == 'Array':
		NumLoci = 0
		LocusGroup = 1
		OutScript = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+str(LocusGroup)+"_mb_genetrees\n#SBATCH -t 48:00:00\n#SBATCH -n 1\nmodule load mrbayes\ncd "+OutFolder+"\n"]
		OutScriptName = OutFolder+OutFilePre+"Grp_"+str(LocusGroup)+"_mb_script.sh"
		for Locus in GoodSeqs:
			NumLoci += 1
			OutFileName = OutFolder+OutFilePre+Locus+"_mb.fa"
			AlignGaps = 0
			for Ind in GoodInds:
				#The version if you want to add blank lines for the missing individuals
				'''
				Record1 = SeqRecord(seq=Seq(OutSeqDict[Locus][Ind]), id = Ind, description = "")
				try:
					AlignGaps.append(Record1)
				except AttributeError:
					AlignGaps = MultipleSeqAlignment([Record1])
				'''
				#The version if you want to leave them out:
				try:
					Record1 = SeqRecord(seq=Seq(SeqIndDict[Locus][Ind]), id = Ind, description = "")
					try:
						AlignGaps.append(Record1)
					except AttributeError:
						AlignGaps = MultipleSeqAlignment([Record1])
				except KeyError:
					"do not add that sequence, because it doesn't exist"
			GapPosDict = GapFinder(Method, AlignGaps, LongSeqList)
			AlignNoGaps = AlignmentGapRemoving(AlignGaps, GapPosDict)
			(AlignOut, SeqsUsing) = MissingIndRemoving(AlignNoGaps)
			if len(AlignNoGaps) != len(AlignOut):
				print("%d sequences that were only gaps were removed." % (len(AlignNoGaps)-len(AlignOut)))
			AlignIO.write(AlignOut, OutFileName, "fasta")
			Line = ScriptFolder+"fasta_to_nexus_mb.py "+OutFilePre+Locus+"_mb.fa\n"
			Line += "mb < "+OutFilePre+Locus+"_mb_mb.nex > "+OutFilePre+Locus+"_mb.log\n"
			OutScript.append(Line)
			if NumLoci == 3:
				OutFileWriting(OutScriptName, OutScript)
				OutLine = "sbatch "+OutScriptName
				SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
				LocusGroup += 1
				NumLoci = 0
				OutScript = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+str(LocusGroup)+"_mb_genetrees\n#SBATCH -t 48:00:00\n#SBATCH -n 1\nmodule load mrbayes\ncd "+OutFolder+"\n"]
				OutScriptName = OutFolder+OutFilePre+"Grp_"+str(LocusGroup)+"_mb_script.sh"
		if NumLoci != 0:
			OutFileWriting(OutScriptName, OutScript)
			OutLine = "sbatch "+OutScriptName
			SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		print("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		sys.stderr.write("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		print("And %d batch scripts were submitted to analyze these loci with MrBayes.\n" % (LocusGroup))
		sys.stderr.write("And %d batch scripts were submitted to analyze these loci with MrBayes.\n" % (LocusGroup))
elif OutMode == 'separate':
	if Mode == 'Parallel':
		OutList1 = [ "#! /bin/bash\n\n"]
		OutList2 = [ ]
		if OGName == "none":
			print("No outgroup will be specified.\n")
			sys.stderr.write("No outgroup will be specified.\n")
		else:
			print("The outgroup is %s.\n" % (OGName))
			sys.stderr.write("The outgroup is %s.\n" % (OGName))
		for Locus in GoodSeqs:
			print Locus
			OutFileName = OutFolder+OutFilePre+Locus+".fa"
			AlignGaps = 0
			LocusIndList = [ ]
			for Ind in GoodInds:
				try:
					Record1 = SeqRecord(seq=Seq(SeqIndDict[Locus][Ind]), id = Ind, description = "")
					try:
						AlignGaps.append(Record1)
					except AttributeError:
						AlignGaps = MultipleSeqAlignment([Record1])
					LocusIndList.append(Ind)
				except KeyError:
					"do not add that sequence, because it doesn't exist"
			GapPosDict = GapFinder(Method, AlignGaps, LongSeqList)
			AlignNoGaps = AlignmentGapRemoving(AlignGaps, GapPosDict)
			(AlignOut, SeqsUsing) = MissingIndRemoving(AlignNoGaps)
			if len(AlignNoGaps) != len(AlignOut):
				print("%d sequences that were only gaps were removed." % (len(AlignNoGaps)-len(AlignOut)))
			AlignIO.write(AlignOut, OutFileName, "fasta")
			Line = "rm "+OutFolder+"RAxML_*."+OutFilePre+Locus+"\n"
			OutList1.append(Line)
			if len(LocusIndList) > 3:
				Line = ScriptFolder+"fasta_to_phylip.py "+OutFileName+" && "
				if OGName == "none":
					Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+".phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -w "+OutFolder+"\n"
				elif OGName == "tree":
					(OGListLocus, LocusOG, NumRootedTree) = OGNamefromTree(SppTree, SeqsUsing)
					if len(SeqsUsing) != NumRootedTree:
						sys.exit("ERROR!!!!!  For locus %s, we started out with %d individuals and %d individuals were in the final tree.  These numbers should be the same!!!!" % (Locus, len(LocusIndList), NumRootedTree))
					Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+".phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+LocusOG+" -w "+OutFolder+"\n"
				else:
					Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+".phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+OGName+" -w "+OutFolder+"\n"
				OutList2.append(Line)
		print("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		sys.stderr.write("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		OutFileName1 = OutFolder+OutFilePre+"analysis_script_raxml1.sh"
		OutFileName2 = OutFolder+OutFilePre+"analysis_script_raxml2.sh"
		Line = "cat "+OutFileName2+" | parallel --joblog "+OutFolder+OutFilePre+"separate_trees.log\n"
		OutList1.append(Line)
		Line = "mkdir "+OutFolder+OutFilePre+"raxmlbs\ncp "+OutFolder+"RAxML_bootstrap."+OutFilePre+"* "+OutFolder+OutFilePre+"raxmlbs\n"
		Line += "cd "+OutFolder+OutFilePre+"raxmlbs\nls RAxML_bootstrap.* > bootstrap_filelist.txt\n"
		Line += "cat "+OutFolder+"RAxML_bestTree."+OutFilePre+"* > "+OutFolder+OutFilePre+"raxmlbs/"+OutFilePre+"bestTrees.tre\n"
		Line += "tar -cf "+OutFolder+OutFilePre+"raxmlbs.tar "+OutFolder+OutFilePre+"raxmlbs\ngzip "+OutFolder+OutFilePre+"raxmlbs.tar.gz\n"
		OutList1.append(Line)
		OutFileWriting(OutFileName1, OutList1)
		OutFileWriting(OutFileName2, OutList2)
	elif Mode == 'Array':
		NumLoci = 0
		LocusGroup = 1
		SBatchList = [ ]
		OutScript = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+str(LocusGroup)+"_genetrees\n#SBATCH -t 10:00:00\n#SBATCH -n 1\nmodule load raxml\n"]
		OutScriptName = OutFolder+OutFilePre+"Grp_"+str(LocusGroup)+"_script.sh"
		for Locus in GoodSeqs:
			print Locus
			OutFileName = OutFolder+OutFilePre+Locus+".fa"
			AlignGaps = 0
			LocusIndList = [ ]
			for Ind in GoodInds:
				try:
					Record1 = SeqRecord(seq=Seq(SeqIndDict[Locus][Ind]), id = Ind, description = "")
					try:
						AlignGaps.append(Record1)
					except AttributeError:
						AlignGaps = MultipleSeqAlignment([Record1])
					LocusIndList.append(Ind)
				except KeyError:
					"do not add that sequence, because it doesn't exist"
			GapPosDict = GapFinder(Method, AlignGaps, LongSeqList)
			AlignNoGaps = AlignmentGapRemoving(AlignGaps, GapPosDict)
			(AlignOut, SeqsUsing) = MissingIndRemoving(AlignNoGaps)
			if len(AlignNoGaps) != len(AlignOut):
				print("%d sequences that were only gaps were removed." % (len(AlignNoGaps)-len(AlignOut)))
			AlignIO.write(AlignOut, OutFileName, "fasta")
			Line = "rm "+OutFolder+"RAxML_*."+OutFilePre+Locus+"\n"
			OutScript.append(Line)
			if len(LocusIndList) > 3:
				NumLoci += 1
				Line = ScriptFolder+"fasta_to_phylip.py "+OutFileName+"\n"
				if OGName == "none":
					Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+".phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -w "+OutFolder+"\n"
				elif OGName == "tree":
					(OGListLocus, LocusOG, NumRootedTree) = OGNamefromTree(SppTree, SeqsUsing)
					if len(SeqsUsing) != NumRootedTree:
						sys.exit("ERROR!!!!!  For locus %s, we started out with %d individuals and %d individuals were in the final tree.  These numbers should be the same!!!!" % (Locus, len(LocusIndList), NumRootedTree))
					Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+".phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+LocusOG+" -w "+OutFolder+"\n"
				else:
					Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+".phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -o "+OGName+" -w "+OutFolder+"\n"
				OutScript.append(Line)
			if NumLoci == 5:
				OutFileWriting(OutScriptName, OutScript)
				OutLine = "sbatch "+OutScriptName
				SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
				JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
				SBatchList.append(JobIDPrev)
				LocusGroup += 1
				NumLoci = 0
				OutScript = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+str(LocusGroup)+"_genetrees\n#SBATCH -t 6:00:00\n#SBATCH -n 1\n\nmodule load raxml\n"]
				OutScriptName = OutFolder+OutFilePre+"Grp_"+str(LocusGroup)+"_script.sh"
		if NumLoci != 0:
			OutFileWriting(OutScriptName, OutScript)
			OutLine = "sbatch "+OutScriptName
			SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
			JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
			SBatchList.append(JobIDPrev)
		print("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		sys.stderr.write("%d sequence files were written, with names such as %s.\n" % (len(GoodSeqs), OutFileName))
		print("And %d batch scripts were submitted to analyze these loci with raxml.\n" % (LocusGroup))
		sys.stderr.write("And %d batch scripts were submitted to analyze these loci with raxml.\n" % (LocusGroup))
		OutScript = ["#! /bin/bash\n#SBATCH -J "+OutFilePre+"combining_genetrees\n#SBATCH -t 1:00:00\n#SBATCH -n 1\n\n"]
		Line = "mkdir "+OutFolder+OutFilePre+"raxmlbs\ncp "+OutFolder+"RAxML_bootstrap."+OutFilePre+"* "+OutFolder+OutFilePre+"raxmlbs\n"
		Line += "cd "+OutFolder+OutFilePre+"raxmlbs\nls RAxML_bootstrap.* > bootstrap_filelist.txt\n"
		Line += "cat "+OutFolder+"RAxML_bestTree."+OutFilePre+"* > "+OutFolder+OutFilePre+"raxmlbs/"+OutFilePre+"bestTrees.tre\n"
		Line += "cd ..\n"
		Line += "tar -cf "+OutFilePre+"raxmlbs.tar "+OutFilePre+"raxmlbs\ngzip "+OutFilePre+"raxmlbs.tar\n"
		OutScript.append(Line)
		OutScriptName = OutFolder+OutFilePre+"combining.sh"
		OutFileWriting(OutScriptName, OutScript)
		OutLine = "sbatch -d afterok:"+":".join(SBatchList)+" "+OutScriptName
		SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
		print("The script %s will be run after these groups of loci have all been analyzed.\n" % (OutScriptName))
		sys.stderr.write("The script %s will be run after these groups of loci have all been analyzed.\n" % (OutScriptName))