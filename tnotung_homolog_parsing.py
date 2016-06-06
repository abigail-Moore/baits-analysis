#! /usr/bin/env python

#tnotung_homolog_parsing.py version 1.0 2 Feb. 2016
#This script parses the trees produced by Notung and determines which sequences are orthologous and which are paralogous.
#They need to be rooted and in Notung format, e.g., RAxML_bipartitions.ppt_nsbest.rearrange.0.ntg

import sys
import re
from collections import defaultdict, OrderedDict
import dendropy
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment


Usage='''
tnotung_homolog_parsing.py version 1.0 2 Feb. 2016
This script parses the homolog table produced by Notung and determines which 
sequences are orthologous and which are paralogous.
tnotung_homolog_parsing.py
[folder in which the Notung trees are found]
[prefix for the Notung trees (anything between RAxML_bipartitions. and the 
locus name)]
[suffix for the Notung trees (anything between the locus name and 
.rearrange.0.ntg]
[folder in which the alignments are found]
[prefix for the alignments (anything before the locus name)]
[suffix for the alignments (anything between the locus name and .fa)]
[file containing the list of loci]
[file containing the group (tab) individual name--to use to say which paralogs 
are present in which/all groups]
[output folder]
[prefix for output files]
[suffix for output files--sequence files only]

For all but the first folder name, substituting "same" causes them to be the
same as the first folder.
For all prefixes and suffixes, substituting "none" causes there to be no 
prefix/suffix.
'''

'''
tnotung_homolog_parsing.py InFolder InFilePre InFilePost AlFolder AlFilePre AlFilePost LocusListFileName IndGroupFileName OutFolder OutFilePre OutFilePost
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 12:
	sys.exit("ERROR!  This script requires 11 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
if InFolder[-1] != "/":
	InFolder += "/"
InFilePre = sys.argv[2]
if InFilePre == "none":
	InFilePre = ""
InFilePost = sys.argv[3]
if InFilePost == "none":
	InFilePost = ""
AlFolder = sys.argv[4]
if AlFolder ==  "same":
	AlFolder = InFolder
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[5]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[6]
if AlFilePost == "none":
	AlFilePost = ""
LocusListFileName = sys.argv[7]
IndGroupFileName = sys.argv[8]
OutFolder = sys.argv[9]
if OutFolder ==  "same":
	OutFolder = InFolder
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[10]
if OutFilePre == "none":
	OutFilePre = ""
OutFilePost = sys.argv[11]
if OutFilePost == "none":
	OutFilePost = ""

ShowTrees = False
Vociferous = False
#MinRename = False

#ShowTrees = True
#Vociferous = True
MinRename = True
#The contigs from minimo have = signs in them, which Dendropy can't read in a Newick format.  So I need to rename the sequences to get rid of those.
#In the next iterations, the contigs will have been renamed from the beginning, so I won't have this problem.
#####################################################################################

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
	print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	#sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is LocusList

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
	#This is IndDict

#GetIndName gets the name of the individual from a sequence, which can either be a contig or a full sequence
#original
def GetIndName(NodeName):
	#for full sequences:
	if (len(str(NodeName).split(".")) > 2):
		IndNameTemp = str(NodeName).split(".")[0]
	#for contigs
	else:
		IndNameTemp = str(NodeName).split("-")[0]
	return IndNameTemp
	#This is IndName

#GetParName gets the name of the old paralog from the sequence name
#The sequences can either be contigs or full sequences
#original
def GetParName(SeqName):
	#for full sequences
	if (len(str(SeqName).split(".")) > 3):
		PName = SeqName.split('.')[1].split('_P')[0]
	#for contigs
	else:
		PName = SeqName.split('_exons_')[-1].split('_P')[0]
	return PName

#SeqLenFinder finds the length of a sequence in an alignment, while ignoring all gaps
#original
def SeqLenFinder(Seq):
	Len = 0
	for Base in Seq:
		if Base != "-":
			Len += 1
	return Len

#NotungFileReading reads the annotated gene tree from a Notung output file and labels the nodes
#according to whether or not a duplication occurred there (Y for duplication, N for no duplication)
#original
def NotungFileReading(FileName):
	try:
		TreeFile = open(FileName, 'rU')
		LineNum = 1
		for Line in TreeFile:
			#The tree is in the first line.
			if LineNum == 1:
				#Various regular expressions to get rid of stuff we don t need and to make the correct node labels.
				r1 = re.compile(r"[rn]\d+:(\d+\.[\dE-]+)\[&&NHX:S=[\w\.]+:Nset=<[\w@\.]+>:H=\w:D=(\w)[:B=0-9\.\]]+")
				r3= re.compile(r"[rn]\d+\[&&NHX:S=[\w\.]+:Nset=<[\w@\.]+>:H=\w:D=(\w)\]")
				r5 = re.compile(r"\[&&NHX:S=[\w\.]+:Nset=<[\w@\.]+>:H=\w+\]")
				TreeString = r5.sub(r"",r3.sub(r"\1",r1.sub(r"\2:\1",Line.rstrip())))
				if MinRename == True:
					r6 = re.compile(r"=")
					TreeString = r6.sub(r"-",TreeString)
			LineNum += 1
		TreeFile.close()
		#turning the string into a tree
		TempTree = dendropy.Tree.get_from_string(TreeString, schema='newick', preserve_underscores=True)
		if ShowTrees == True: print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	except IOError:
		TempTree = "nothing"
	return TempTree
	#This is LocusTree_orig

#SeqFileReading reads a sequence file and puts the sequences in a dictionary.
#from tparcomb_final.py (although it was maybe from somewhere else before that)
def SeqFileReading(FileName, SeqFormat):
	DictTemp = { }#DictTemp[SeqName] = Seq
	InFile = open(FileName, 'rU')
	for record in SeqIO.parse(InFile, SeqFormat):
		if MinRename == True:
			MinRe = re.compile("=")
			SeqName = MinRe.sub("-",record.id)
			DictTemp[SeqName] = str(record.seq)
		else:
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

#DuplicationFinder goes through each node on the the Notung tree and determines whether an actual duplication has occurred
#original
def DuplicationFinder(TempTree):
	if Vociferous == True: print("Resolving duplications due to incongruences between the gene tree and the species tree.\n")
	NDs = 0
	NDNs = 0
	for Node in TempTree.preorder_node_iter():
		if Node.label == "Y":
			NDs += 1
			ParNum = 1
			DictTemp = defaultdict(list)
			#if it did, finding out which individuals are daughters of that node
			for Child_Node in Node.child_node_iter():
				for Daughter_Node in Child_Node.leaf_nodes():
					DictTemp[ParNum].append(GetIndName(Daughter_Node.taxon.label))
				ParNum += 1
			#if none of the individuals are present twice, there probably isn't really a duplication at that node.
			#It is more likely that the species and gene trees are incongruent.
			OVSet = set(DictTemp[1]).intersection(DictTemp[2])
			if len(OVSet) == 0:
				if Vociferous == True: print("All taxa present only once!")
				#so changing the node label so that there is no duplication
				Node.label = "N"
			#numbering the duplication, if we accept it
			else:
				NDNs += 1
				Node.label = "Y"+str(NDNs)
		elif not Node.label:
			Node.label = "N"
		elif Node.label != "N":
			sys.exit("ERROR!!  This tree was incorrectly parsed as the node label %s was found.\n" % (Node.label))
	if ShowTrees == True: print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	if Vociferous == True: print("Prior to resolving duplications caused by differences between the gene tree and the species tree, there were %d duplications.  Now there are %d.\n" % (NDs, NDNs))
	return TempTree
	#These are LocusTree_new, NumDups, NumDupsNew.

#SeqCombiner goes through the tree and examines each putative duplication in more detail to determine whether it is really a duplication
#or whether some of the sequences can be merged.
#original
def SeqCombiner(TempTree, InSeqDict):
	if Vociferous == True: print("Resolving duplications that have sequences that should have been combined into one sequence.\n")
	NumDupsOld = 0
	NumDupsNew = 0
	for Node in TempTree.postorder_node_iter():
		if Node.label[0] == "Y":
			NumDupsOld += 1
			NumDupsNew += 1
			DaughterDict = defaultdict(list)
			DaughterList = [ ]
			RemoveList = [ ]
			RenameDict = { }
			CombinableInds = [ ]
			#determine which individuals are present more than once
			for Daughter_Node in Node.leaf_nodes():
				IndName = GetIndName(Daughter_Node.taxon.label)
				DaughterDict[IndName].append(Daughter_Node.taxon.label)
				DaughterList.append(IndName)
			DaughterList = list(set(DaughterList))
			DupInds = [ ]
			for IndName in DaughterDict:
				if len(DaughterDict[IndName]) > 1:
					DupInds.append(IndName)
					#if Vociferous == True: print("Individual %s is present more than once at node %s.\n" % (IndName, Node.label))
			if Vociferous == True: print("Of the %d individuals at node %s, the following %d were present more than once: %s.\n" % (len(DaughterList), Node.label, len(DupInds), ", ".join(DupInds)))
			#Check first to determine whether there are fewer than xx individuals or less than xx percent of the total?
			if (len(DupInds) < 5) or (float(len(DupInds))/len(DaughterList) < 0.4):
				if Vociferous == True: print("There were few enough duplicated individuals, that we will try to combine their sequences.\n")
				#for each of those individuals, try to combine its sequences
				for IndName in DupInds:
					#if Vociferous == True: print("Individual %s has the following sequences: %s.\n" % (IndName, ", ".join(DaughterDict[IndName])))
					(ConSeq, NumAmbig, AmbigSitesList, NumSeqs, Overlap) = ConSeqMaker(InSeqDict, DaughterDict[IndName])
					if NumAmbig < 5:
						CombinableInds.append(IndName)
						ParNameList =sorted(list(set([GetParName(SeqName) for SeqName in DaughterDict[IndName]])))
						#Getting the NewParName, and trying to avoid Ambig, if possible
						NewParName = ""
						for ParName in ParNameList:
							if ParName[:5] != "Ambig":
								NewParName = ParName
						#but, if all of the potential paralogs are Ambig, then choosing the first one
						if NewParName == "":
							NewParName = ParNameList[0]
						NumSeqsTot = NumSeqs
						NumAmbigsTot = NumAmbig
						for SeqName in DaughterDict[IndName]:
							if len(SeqName.split("seqs")) == 2:
								NumSeqsTot += int(SeqName.split("seqs")[0].split(".")[-1])-1
								NumAmbigsTot += int(SeqName.split("ambig")[0].split(".")[-1])
						NewSeqName = ".".join([IndName, ParNameList[0], str(NumSeqsTot)+"seqs", str(NumAmbigsTot)+"ambig", "len"+str(SeqLenFinder(ConSeq))])
						InSeqDict[NewSeqName] = ConSeq
						#deleting the shorter sequences from the tree
						LengthLongest = 0
						LongestSeq = ""
						for SeqName in DaughterDict[IndName]:
							SeqLen = SeqLenFinder(InSeqDict[SeqName])
							if SeqLen > LengthLongest:
								LengthLongest = SeqLen
								LongestSeq = SeqName
						for SeqName in DaughterDict[IndName]:
							#remove shorter sequences
							if SeqName != LongestSeq:
								RemoveList.append(SeqName)
							#rename longer sequence
							else:
								RenameDict[SeqName] = NewSeqName
							#in either case, remove it from the InSeqDict
							del InSeqDict[SeqName]
						if Vociferous == True: print("The sequences (%s) for individual %s were combinable into %s!\n" % (", ".join(DaughterDict[IndName]), IndName, NewSeqName))
					else:
						if Vociferous == True: print("The sequences for individual %s were not combinable.  :(\n" % (IndName))
			#going through the leaf nodes again and removing/renaming them, if necessary
			for Daughter_Node in Node.leaf_nodes():
				OldName = Daughter_Node.taxon.label
				if OldName in RenameDict.keys():
					Daughter_Node.taxon.label = RenameDict[OldName]
					#if Vociferous == True: print("Taxon %s renamed to %s.\n" % (OldName, RenameDict[OldName]))
				elif OldName in RemoveList:
					TempTree.prune_taxa([Daughter_Node.taxon], update_bipartitions=False)
					#if Vociferous == True: print("Sequence %s removed.\n" % (Daughter_Node.taxon.label))
			#relabeling the node
			if sorted(CombinableInds) == sorted(DupInds):
				if Vociferous == True: print("Duplication %s successfully resolved!\n" % (Node.label))
				Node.label = "N"
				NumDupsNew = NumDupsNew -1
	if ShowTrees == True: print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	if Vociferous == True: print("Prior to combining separated parts of sequences of the same paralog, there were %d duplications.  Now there are %d.\n" % (NumDupsOld, NumDupsNew))
	return (TempTree, InSeqDict)

#DuplicationMapper goes through the tree and figures out the nodes above and below each duplication, for subsequent mapping\
#original
def DuplicationMapper(TempTree):
	DictTemp = defaultdict(dict)#DictTemp[Node.label]['sister'/'daughters'] = list of individuals 
	for Node in TempTree.preorder_node_iter():
		#look at each duplication
		if Node.label[0] == "Y":
			DictTemp[Node.label]['daughters'] = list(set([GetIndName(Daughter_Node.taxon.label) for Daughter_Node in Node.leaf_nodes()]))
			DictTemp[Node.label]['sister'] = [ ]
			for Sis_Node in Node.sister_nodes():
				DictTemp[Node.label]['sister'] += list(set([GetIndName(Niece_Node.taxon.label) for Niece_Node in Sis_Node.leaf_nodes()]))
	return DictTemp
	#This is DupPosDict[Locus]

#IndDupsRemoving goes through the tree, finds duplications that consist of only one individual, and keeps only the longest
#sequence in those duplications
#original
def IndDupsRemoving(TempTree, InSeqDict):
	if Vociferous == True: print("Resolving duplications that consist of a single individual.\n")
	NumDupsOld = 0
	NumDupsNew = 0
	for Node in TempTree.postorder_node_iter():
		#look at each duplication
		if Node.label[0] == "Y":
			NumDupsOld += 1
			NumDupsNew += 1
			#figure out which individuals are present
			DaughterList = [Daughter_Node.taxon.label for Daughter_Node in Node.leaf_nodes()]
			IndList = [GetIndName(SeqName) for SeqName in DaughterList]
			IndList = list(set(IndList))
			#if the duplication is within a single individual
			if len(IndList) == 1:
				if Vociferous == True: print("Duplication %s consists only of sequences from individual %s.\n" % (Node.label, IndList[0]))
				#find the longest sequence
				LengthLongest = 0
				LongestSeq = ""
				for SeqName in DaughterList:
					SeqLen = SeqLenFinder(InSeqDict[SeqName])
					if SeqLen > LengthLongest:
						LengthLongest = SeqLen
						LongestSeq = SeqName
				#remove the remaining sequence(s)
				for Daughter_Node in Node.leaf_nodes():
					if Daughter_Node.taxon.label != LongestSeq:
						TempTree.prune_taxa([Daughter_Node.taxon], update_bipartitions=False)
						if Vociferous == True: print("For node %s, sequence %s will be removed and sequence %s will be kept.\n" % (Node.label, Daughter_Node.taxon.label, LongestSeq))
				Node.label = "N"
				NumDupsNew = NumDupsNew -1
			else:
				if Vociferous == True: print("Duplication %s consists of sequences from %d individuals.\n" % (Node.label, len(IndList)))
	if ShowTrees == True: print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	if Vociferous == True: print("Prior to resolving duplications consisting of a single individual, there were %d duplications.  Now there are %d.\n" % (NumDupsOld, NumDupsNew))
	return (TempTree)			

#ParalogNamer gives a consistent name to each of the new paralogs
#This really assumes that the tree will be bifurcating at the nodes where there is a duplication.  Some of it may work otherwise,
#but I don't think it will work properly.
#original
def ParalogNamer(TempTree, Locus):
	#going through the tree with the corrected node labels and splitting it at duplications:
	DictTemp = defaultdict(dict)#DictTemp[DupNum][ParNum]['NestedPars']['Seqs']['ParNames'] = lists of things
	OutDictTemp = { }#OutDictTemp[SeqName] = NewParName
	MaxParDict = { }#MaxParDict[ParName] = highest integer it has been given
	RDict = { }#RDict[ParName] = ParRank
	RDict[Locus] = 0
	ParRank = 1
	AmbigParNum = 1
	for Node in TempTree.preorder_node_iter():
		if Node.label[0] == "Y":
			DupNum = Node.label[1:]
			ParNum = 0
			DupNamesList = [ ]
			#for each half of the duplication,
			for Child_Node in Node.child_node_iter():
				ParNum += 1
				#making the dictionary of things to fill out for this half
				DictTemp[DupNum][ParNum] = defaultdict(list)
				#looking for nested duplications
				if Child_Node.label[0] == "Y":
					DictTemp[DupNum][ParNum]['NestedPars'].append(Child_Node.label[1:])
				for Daughter_Node2 in Child_Node.preorder_iter():
					if Daughter_Node2.label[0] == "Y":
						DictTemp[DupNum][ParNum]['NestedPars'].append(Daughter_Node2.label[1:])
				for Daughter_Node in Child_Node.leaf_nodes():
					#making a list of sequences/OTUs in that half
					SeqName = Daughter_Node.taxon.label
					DictTemp[DupNum][ParNum]['Seqs'].append(SeqName)
					#and a list of the paralog names those OTUs already have
					#they can potentially have a name because they were given one from an earlier paralog
					try:
						ParNameTemp = OutDictTemp[SeqName]
						if ParNameTemp[0] != "?":
							DictTemp[DupNum][ParNum]['ParNames'].append(OutDictTemp[SeqName])
						else:
							DictTemp[DupNum][ParNum]['ParNames'].append(GetParName(SeqName))
					#or they can still have their original name
					except KeyError:
						DictTemp[DupNum][ParNum]['ParNames'].append(GetParName(SeqName))
				#condensing the list of paralog names for that half of the duplication
				ListTemp = DictTemp[DupNum][ParNum]['ParNames']
				ListTemp = list(set(ListTemp))
				DictTemp[DupNum][ParNum]['ParNames'] = ListTemp
				DupNamesList += ListTemp
				if Vociferous == True: print("Half %d of duplication %s had %d nested duplications and %d sequences, and already had daughter paralogs with the following names: %s\n" % (ParNum, DupNum, len(DictTemp[DupNum][ParNum]['NestedPars']), len(DictTemp[DupNum][ParNum]['Seqs']), ", ".join(DictTemp[DupNum][ParNum]['ParNames'])))
			#condensing the set of names for the duplication
			DupNamesList = list(set(DupNamesList))
			#this is the final number, so I can iterate over ParNum below
			NumPars = ParNum
			if Vociferous == True: print("Duplication number %s has %d nested duplications and already has daughter paralogs with the following names: %s.\n" % (DupNum, len(DictTemp[DupNum][1]['NestedPars'])+len(DictTemp[DupNum][2]['NestedPars']), ", ".join(DupNamesList))) 
			#if each half of the duplication had a single name:
			if (len(DictTemp[DupNum][1]['ParNames']) == 1) and (len(DictTemp[DupNum][2]['ParNames']) == 1):
				#the names can be different
				if len(DupNamesList) == ParNum:
					for ParNum in DictTemp[DupNum]:
						NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]
						DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
						#if that name already exists,then it needs to get that name with a number
						if NewParName in RDict.keys():
							try:
								NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]+"_"+str(MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]]+1)
								MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]] += 1
							except KeyError:
								NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]+"_1"
								MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]] = 1
						#if not, it can just get that name
						RDict[NewParName] = ParRank
						ParRank += 1
						for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
							OutDictTemp[SeqName] = NewParName
				#or the names can be the same
				else:
					#We need to go through this in order for the assigning of paralog numbers to work
					for ParNum in sorted(DictTemp[DupNum].keys()):
						try:
							NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]+"_"+str(MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]]+1)
							MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]] += 1
						except KeyError:
							NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]+"_1"
							MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]] = 1
						DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
						RDict[NewParName] = ParRank
						ParRank += 1
						for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
							OutDictTemp[SeqName] = NewParName
						if Vociferous == True: print("The new paralog name for half %d of duplication number %s is %s.\n" % (ParNum, DupNum, ", ".join(DictTemp[DupNum][ParNum]['NewParNames'])))
			#It is more complicated if this is not the case.
			else:
				if Vociferous == True: print("The case is more complicated.")
				for ParNum in DictTemp[DupNum]:
					#maybe one of the halves is alright
					#if that half has one name,
					if len(DictTemp[DupNum][ParNum]['ParNames']) == 1:
						ParNameTemp = DictTemp[DupNum][ParNum]['ParNames'][0]
						#the name can only belong to that half, in which case it can be the name of the paralog
						if ((ParNum == 1) and (ParNameTemp not in DictTemp[DupNum][2]['ParNames'])) or ((ParNum == 2) and (ParNameTemp not in DictTemp[DupNum][1]['ParNames'])):
							if Vociferous == True: print("Half %d of duplication %s has a unique paralog name: %s.\n" % (ParNum, DupNum, ParNameTemp))
							NewParName = ParNameTemp
							#if that name already exists,then it needs to get that name with a number
							if NewParName in RDict.keys():
								try:
									NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]+"_"+str(MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]]+1)
									MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]] += 1
								except KeyError:
									NewParName = DictTemp[DupNum][ParNum]['ParNames'][0]+"_1"
									MaxParDict[DictTemp[DupNum][ParNum]['ParNames'][0]] = 1
							#if not, it can just get that name
							RDict[NewParName] = ParRank
							ParRank += 1
							DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
							for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
								OutDictTemp[SeqName] = NewParName
						#or else that name can also be present in the other half, in which case it becomes the paralog name but with an additional integer.
						else:
							try:
								NewParName = ParNameTemp+"_"+str(MaxParDict[ParNameTemp]+1)
								MaxParDict[ParNameTemp] += 1
							except KeyError:
								NewParName = ParNameTemp+"_1"
								MaxParDict[ParNameTemp] = 1
							DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
							RDict[NewParName] = ParRank
							ParRank += 1
							for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
								OutDictTemp[SeqName] = NewParName
							if Vociferous == True: print("The new paralog name for half %d of duplication number %s is %s.\n" % (ParNum, DupNum, NewParName))
					else:
						#If we have multiple paralog names, I want to rename the sequences, but flag them as uncertain.
						#Hopefully this will be resolved as we go up the tree.
						#If not, then we can go back and correct it later.
						try:
							NewParName = "?"+Locus+"_"+str(MaxParDict[Locus]+1)
							MaxParDict[Locus] += 1
						except KeyError:
							NewParName = Locus+"_1"
							MaxParDict[Locus] = 1
						DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
						for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
							OutDictTemp[SeqName] = NewParName
						RDict[NewParName] = ParRank
						ParRank += 1
						if Vociferous == True: print("The new paralog name for half %d of duplication number %s is uncertain, but is temporarily %s.\n" % (ParNum, DupNum, NewParName))
	#Then going through all of the SeqNames in the tree to make sure they have an associated paralog in the dictionary
	SeqsNotClassif = 0
	#go through the duplications, to see if their sequences were properly classified
	for DupNum in DictTemp:
		for ParNum in DictTemp[DupNum]:
			NewParList = DictTemp[DupNum][ParNum]['NewParNames']
			for ParName in NewParList:
				if ParName[0] == '?':
					AmbigParSeqs = [ ]
					#if ambiguously classified sequences are found,
					for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
						if OutDictTemp[SeqName] == ParName:
							AmbigParSeqs.append(SeqName)
					#remove the ? from their name and update the various lists
					DictTemp[DupNum][ParNum]['NewParNames'].remove(ParName)
					if len(AmbigParSeqs) != 0:
						NewParName = ParName[1:]
						DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
						for SeqName in AmbigParSeqs:
							if Vociferous == True: print("ERROR!  Sequence %s was not properly classified according to paralog!\n" % (SeqName))
							SeqsNotClassif += 1
							OutDictTemp[SeqName] = NewParName
						NewParRank = RDict[ParName]
						RDict[NewParName] = NewParRank
					del RDict[ParName]
	#Then go through all of the sequences in the tree, to make sure that they have paralogs.
	for LeafNode in TempTree.leaf_nodes():
		SeqName = LeafNode.taxon.label
		#see if the sequence is already has a paralog name
		try:
			ParName = OutDictTemp[SeqName]
		#If not, they were not duplicated, and their paralog should be the locus name.
		except KeyError:
			OutDictTemp[SeqName] = Locus
	print("%s of %s sequences could not be classified correctly.\n" % (SeqsNotClassif, len(OutDictTemp)))
	return DictTemp, OutDictTemp, RDict
	#These are ParInfoDict[Locus], NewParDict[Locus], and RankDict.

#OrthoParFinder goes through a tree and separates it into sets of orthologous genes, while making separate groupings for the
#various paralogs
#original
def OrthoParFinder(TempTree, PNDict, RDict, LName):
	DictTemp = defaultdict(dict)#DictTempSeqs[Node.label]['SeqNames'/'ParList'/'Name']
	for Node in TempTree.postorder_node_iter():
		if Node.label[0] == "Y":
			NodeDict = defaultdict(list)
			NodeSeqDict = defaultdict(list)
			for Child_Node in Node.child_node_iter():
				for Daughter_Node in Child_Node.leaf_nodes():
					#making a list of individuals in that half
					SeqName = Daughter_Node.taxon.label
					NodeDict[Child_Node].append(GetIndName(SeqName))
					NodeSeqDict[Child_Node].append(SeqName)
			MostInds = 0
			BestNode = ""
			#choose the node with the most individuals (or the first one, if they both have the same number of individuals)
			for Child_Node in NodeDict:
				NumInds = len(list(set(NodeDict[Child_Node])))
				if NumInds > MostInds:
					MostInds = NumInds
					BestNode = Child_Node
			#cut the tree on the branch leading to the other node
			for Child_Node in NodeDict:
				if Child_Node != BestNode:
					if Vociferous == True: print("For node %s, the subtree with the following sequences will be pruned %s.\n" % (Node.label, ", ".join(NodeSeqDict[Child_Node])))
					DictTemp[Node.label]['SeqNames'] = NodeSeqDict[Child_Node]
					ParList = list(set([PNDict[Seq] for Seq in NodeSeqDict[Child_Node]]))
					DictTemp[Node.label]['ParList'] = ParList
					if Vociferous == True: print("This subtree includes the following paralogs: %s.\n" % (", ".join(ParList)))
					LowestRank = 5000
					BestPar = ""
					for ParName in ParList:
						Rank = RDict[ParName]
						if (Rank < LowestRank) and (Rank != 0):
							LowestRank = Rank
							BestPar = ParName
					if Vociferous == True: print("The new name for this paralog will be %s.\n" % (BestPar))
					DictTemp[Node.label]['Name'] = BestPar
					TempTree.prune_subtree(Child_Node, update_bipartitions=False)
					#print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	DictTemp['root']['SeqNames'] = [LeafNode.taxon.label for LeafNode in TempTree.leaf_nodes()]
	ParList = list(set([(PNDict[Seq]) for Seq in DictTemp['root']['SeqNames']]))
	DictTemp['root']['ParList'] = ParList
	if Vociferous == True: print("The largest group of sequences includes the following %d individuals: %s\n" % (len(DictTemp['root']['SeqNames']), ", ".join(DictTemp['root']['SeqNames'])))
	if Vociferous == True: print("It includes the following paralogs: %s.\n" % (", ".join(ParList)))
	LowestRank = 5000
	BestPar = ""
	for ParName in ParList:
		Rank = RDict[ParName]
		if (Rank < LowestRank) and (Rank != 0):
			LowestRank = Rank
			BestPar = ParName
	if BestPar == "":
		BestPar = LName
	if Vociferous == True: print("The new name for this paralog will be %s.\n" % (BestPar))
	DictTemp['root']['Name'] = BestPar
	return DictTemp

#SeqsPerIndFinder looks through a tree and makes a dictionary with the number of sequences that belong to each individual
#original
def SeqsPerIndFinder(TempTree):
	DictTemp = defaultdict(int)
	for LeafNode in TempTree.leaf_node_iter():
		IndName = GetIndName(LeafNode.taxon.label)
		DictTemp[IndName] += 1
	return DictTemp
	#These are the portions of SeqsPerInd and SeqsPerIndPP for each locus. 

#ConSeqMaker makes a consensus sequence from a group of aligned sequences
#from tcontigs_to_fixed_paralogs.py
def ConSeqMaker(SeqDict,SeqList):
	AlignTemp = [ ]
	#find the sequences for that individual
	for Contig in SeqList:
		SeqTemp = SeqRecord(seq=(Seq(SeqDict[Contig])), id=Contig),
		AlignTemp += SeqTemp
	#put them in an alignment
	AlignTemp = MultipleSeqAlignment(AlignTemp)
	NSeqs = len(AlignTemp)
	#make the consensus of the alignment
	ConSeq = ""
	AmbigNucs = 0
	Overlap = 0
	AmbigNucList = [ ]
	for SeqPos in range(0, len(SeqDict[Contig])):
		PosNucs = [ ]
		NumNucs = 0
		for record in AlignTemp:
			if record[SeqPos] != '-':
				PosNucs += record[SeqPos]
				NumNucs += 1
		if NumNucs > 1:
			Overlap += 1
		PosNucs = list(set(PosNucs))
		if len(PosNucs) == 1:
			ConSeq += PosNucs[0]
		elif len(PosNucs) > 2:
			ConSeq += 'n'
			AmbigNucs += 1
			AmbigNucList.append(SeqPos)
		elif len(PosNucs) == 2:
			if 'n' in PosNucs: ConSeq += 'n'
			elif 'm' in PosNucs:
				if 'a' in PosNucs: ConSeq += 'm'
				elif 'c' in PosNucs: ConSeq += 'm'
				else: ConSeq += 'n'
			elif 'r' in PosNucs:
				if 'a' in PosNucs: ConSeq += 'r'
				elif 'g' in PosNucs: ConSeq += 'r'
				else: ConSeq += 'n'
			elif 'w' in PosNucs:
				if 'a' in PosNucs: ConSeq += 'w'
				elif 't' in PosNucs: ConSeq += 'w'
				else: ConSeq += 'n'
			elif 's' in PosNucs:
				if 'c' in PosNucs: ConSeq += 's'
				elif 'g' in PosNucs: ConSeq += 's'
				else: ConSeq += 'n'
			elif 'y' in PosNucs:
				if 'c' in PosNucs: ConSeq += 'y'
				elif 't' in PosNucs: ConSeq += 'y'
				else: ConSeq += 'n'
			elif 'k' in PosNucs:
				if 'g' in PosNucs: ConSeq += 'k'
				elif 't' in PosNucs: ConSeq += 'k'
				else: ConSeq += 'n' 
			elif 'a' in PosNucs:
				if 'c' in PosNucs: ConSeq += 'm'
				elif 'g' in PosNucs: ConSeq += 'r'
				elif 't' in PosNucs: ConSeq += 'w'
				else:
					print("ERROR!!! We have strange bases in our sequence: %s\n" % (" ".join(PosNucs)))
					ConSeq += 'n'
			elif 'c' in PosNucs:
				if 'g' in PosNucs: ConSeq += 's'
				elif 't' in PosNucs: ConSeq += 'y' 
				else:
					print("ERROR!!! We have strange bases in our sequence: %s\n" % (" ".join(PosNucs)))
					ConSeq += 'n'
			elif 'g' in PosNucs:
				if 't' in PosNucs: ConSeq += 'k'
				else:
					print("ERROR!!! We have strange bases in our sequence: %s\n" % (" ".join(PosNucs)))
					ConSeq += 'n'
			AmbigNucs += 1
			AmbigNucList.append(SeqPos)
		elif len(PosNucs) == 0:
			ConSeq += '-'
	return (ConSeq, AmbigNucs, AmbigNucList, NSeqs, Overlap)

#######################################################################################

LocusList = CaptureColumn(LocusListFileName, 0)
IndGroupDict = DictFromFile(IndGroupFileName, 1, 0)

#The dictionaries that will be filled out:
TreeDict = { }#TreeDict[Locus] = LocusTree--nodes labeled according to duplications, with spurious duplications removed
ParInfoDict = defaultdict(dict)#ParInfoDict[Locus][DupNum][ParNum][various keys] = various lists
NewParDict = defaultdict(dict)#NewParDict[Locus][SeqName] = ParName
NewParInfoDict = defaultdict(dict)#NewParInfoDict[Locus][DupName]['SeqNames'/'ParList'/'Name'] = list of SeqNames/list of ParNames/final name of paralog
DupPosDict = defaultdict(dict)#DupPosDict[Locus][DupName]['sister'/'daughters'] = list of individuals
SeqsPerIndGF = defaultdict(dict)#SeqsPerIndGF[Locus][Ind] = NumSeqs
SeqsPerInd = defaultdict(int)#SeqsPerInd[Ind] = NumSeqs
SeqLenParDict = defaultdict(dict)#SeqLenParDict[Locus][Paralog][Ind] = SeqLen
IndsPerGroupDict = defaultdict(dict)#IndsPerGroupDict[Locus][ParName][GroupName] = NumSeqs per GroupName
#Also dictionaries that include the duplications at the tips, for looking at polyploidy:
DupPosDictPP = defaultdict(dict)#DupPosDictPP[Locus][DupName]['sister'/'daughters'] = list of individuals
SeqsPerIndGFPP = defaultdict(dict)#SeqsPerIndGFPP[Locus][Ind] = NumSeqs

LocusListOut = [ ]
LociRead = 0
LociUnRead = 0
for Locus in LocusList:
	sys.stderr.write(Locus)
	#reading the tree from Notung (.ntg format)
	TreeFileName = InFolder+"RAxML_bipartitions."+InFilePre+Locus+InFilePost+".rearrange.0.ntg"
	LocusTree_orig = NotungFileReading(TreeFileName)
	if LocusTree_orig == "nothing":
		print("No tree file for locus %s.\n" % (Locus))
		sys.stderr.write("No tree file for locus %s.\n" % (Locus))
		LociUnRead += 1
	else:
		LociRead += 1
		#setting up the dictionaries
		SeqsPerIndGF[Locus] = defaultdict(int)
		SeqLenParDict[Locus] = defaultdict(dict)
		#reading the alignment (fasta format)
		AlFileName = AlFolder+AlFilePre+Locus+AlFilePost+".fa"
		LocusAl_orig = SeqFileReading(AlFileName, 'fasta')
		#going through each node in the tree and determining whether Notung reconstructed a duplication at that node
		(LocusTree_new) = DuplicationFinder(LocusTree_orig)
		#Go through the duplications and try to combine the sequences that can be combined
		(LocusTree_pruned, NewSeqDict) = SeqCombiner(LocusTree_new, LocusAl_orig)
		#choose the longer sequence for each duplication involving only one individual that could not be successfully resolved
		(LocusTree_nodups) = IndDupsRemoving(LocusTree_pruned, LocusAl_orig)
		#then figuring out where all the duplications are (first including duplications in the tips)
		DupPosDictPP[Locus] = DuplicationMapper(LocusTree_pruned)
		#then figuring out where all the duplications are (excluding duplications in the tips)
		DupPosDict[Locus] = DuplicationMapper(LocusTree_nodups)
		#finding out how many sequences each individual has (first including sequences that are just duplicated in the tips)
		SeqsPerIndGFPP[Locus] = SeqsPerIndFinder(LocusTree_pruned)
		#finding out how many sequences each individual has (then excluding duplications in the tips)
		SeqsPerIndGF[Locus] = SeqsPerIndFinder(LocusTree_nodups)
		#finding a name for each duplication
		(ParInfoDict[Locus], NewParDict[Locus], RankDict) = ParalogNamer(LocusTree_nodups, Locus)
		#go through and find more or less orthologous sequences, using the half of each duplication with the largest taxon sampling
		NewParInfoDict[Locus] = OrthoParFinder(LocusTree_nodups, NewParDict[Locus], RankDict, Locus)
		#write a new sequence file for the gene family
		OutFileNameGF = OutFolder+OutFilePre+Locus+OutFilePost+".fa"
		#write new sequence files for each new paralog
		for DupNode in NewParInfoDict[Locus]:
			#one copy with the original sequence names
			ParName = NewParInfoDict[Locus][DupNode]['Name']
			LocusListOut.append(Locus+"_"+ParName+"\n")
			OutFileName = OutFolder+OutFilePre+Locus+"_"+ParName+OutFilePost+"_seqnames.fa"
			#making a dictionary with just the sequences for this paralog
			SeqDict = { }
			for SeqName in NewParInfoDict[Locus][DupNode]['SeqNames']:
				SeqDict[SeqName] = NewSeqDict[SeqName]
			SeqFileWriting(OutFileName, SeqDict, 'fasta')
			#and one copy with the sequence just named after the individual, for further analysis
			SeqDictInds = { }
			OutFileName = OutFolder+OutFilePre+ParName+OutFilePost+".fa"
			if ParName == "singlecopy":
				ParName = Locus
			IndsPerGroupDict[Locus][ParName] = defaultdict(int)
			for SeqName in SeqDict:
				IndName = GetIndName(SeqName)
				GroupName = IndGroupDict[IndName]
				SeqDictInds[IndName] = SeqDict[SeqName]
				#also filling out some dictionaries
				SeqsPerInd[IndName] += 1
				SeqLenParDict[Locus][ParName][IndName] =  SeqLenFinder(SeqDict[SeqName])
				IndsPerGroupDict[Locus][ParName][GroupName] += 1
				IndsPerGroupDict[Locus][ParName]['all'] += 1
			SeqFileWriting(OutFileName, SeqDictInds, 'fasta')
		SeqFileWriting(OutFileNameGF, NewSeqDict, 'fasta')

print("Of the %d loci, %d had tree files and %d did not.\n" % (len(LocusList), LociRead, LociUnRead))
sys.stderr.write("Of the %d loci, %d had tree files and %d did not.\n" % (len(LocusList), LociRead, LociUnRead))
##################################################################################
#writing the output files

#making a new pardict for future analyses
OutList = [ ]
for Locus in NewParDict:
	for SeqName in NewParDict[Locus]:
		Line = Locus+"\t"+NewParDict[Locus][SeqName]+"\t"+SeqName+"\t"+SeqName.split("_")[0]+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"pardict.txt"
OutFileWriting(OutFileName, OutList)

#making  list of nested paralogs
OutList = ['Locus\tParalog_Name\tNested_Paralogs\n']
for Locus in NewParInfoDict:
	for DupName in NewParInfoDict[Locus]:
		Line = ("\t").join([Locus, NewParInfoDict[Locus][DupName]['Name'], ",".join(NewParInfoDict[Locus][DupName]['ParList'])])+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"nested_paralog_dict.txt"
OutFileWriting(OutFileName, OutList)

#making a list of duplication positions for eventual mapping
OutList = ['Locus\tDuplication_Name\tDuplicated_Individuals\tSister_Group\n']
for Locus in DupPosDict:
	for DupName in DupPosDict[Locus]:
		Line = ("\t").join([Locus, DupName, ",".join(DupPosDict[Locus][DupName]['daughters']), ",".join(DupPosDict[Locus][DupName]['sister'])])+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"dup_pos_dict.txt"
OutFileWriting(OutFileName, OutList)

#making a list of duplication positions for eventual mapping, and the same including tip duplications
OutList = ['Locus\tDuplication_Name\tDuplicated_Individuals\tSister_Group\n']
for Locus in DupPosDictPP:
	for DupName in DupPosDictPP[Locus]:
		Line = ("\t").join([Locus, DupName, ",".join(DupPosDictPP[Locus][DupName]['daughters']), ",".join(DupPosDictPP[Locus][DupName]['sister'])])+"\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"dup_pos_dict_plus_tips.txt"
OutFileWriting(OutFileName, OutList)

#making a list of loci that had trees
OutFileName = OutFolder+OutFilePre+"locus_list_out.txt"
OutFileWriting(OutFileName, LocusListOut)

IndList = sorted(SeqsPerInd.keys())
LocusKeyList = sorted(SeqsPerIndGF.keys())

#printing the number of sequences per individual: SeqsPerInd
OutFileName = OutFolder+OutFilePre+"seqs_per_ind.txt"
OutList = ['Individual_Name\tTotal_Sequences\n']
for IndName in IndList:
	Line = IndName+"\t"+str(SeqsPerInd[IndName])+"\n"
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#printing the number of sequences per individual per locus: SeqsPerIndGF
OutFileName = OutFolder+OutFilePre+"seqs_per_ind_locus.txt"
OutList = ['Individual_Name'+"\t"+"\t".join(LocusKeyList)+"\n"]
for IndName in IndList:
	Line = IndName
	for Locus in LocusKeyList:
		Line += "\t"+str(SeqsPerIndGF[Locus][IndName])
	Line += "\n"
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)


#printing the number of sequences per individual per locus, including duplications at the tips: SeqsPerIndGFPP
OutFileName = OutFolder+OutFilePre+"seqs_per_ind_locus_plus_tips.txt"
OutList = ['Individual_Name'+"\t"+"\t".join(LocusKeyList)+"\n"]
for IndName in IndList:
	Line = IndName
	for Locus in LocusKeyList:
		Line += "\t"+str(SeqsPerIndGFPP[Locus][IndName])
	Line += "\n"
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#printing a presence/absence table and a table of the length of the sequences that are present for each paralog and individual
OutFileName1 = OutFolder+OutFilePre+"seq_presence_absence.txt"
OutFileName2 = OutFolder+OutFilePre+"seq_length.txt"
OutDict1 = defaultdict(str)
OutDict2 = defaultdict(str)
KeyLine = "IndividualName\t"
for Locus in LocusKeyList:
	PKeyList = sorted(SeqLenParDict[Locus].keys())
	for Paralog in PKeyList:
		for Ind in IndList:
			try:
				OutDict2[Ind] += str(SeqLenParDict[Locus][Paralog][Ind])+"\t"
				OutDict1[Ind] += 'x\t'
			except KeyError:
				OutDict2[Ind] += '\t'
				OutDict1[Ind] += '\t'
	KeyLine += "\t".join(PKeyList)+"\t"
OutList1 = [KeyLine[:-1]+"\n"]
OutList2 = [KeyLine[:-1]+"\n"]
for Ind in IndList:
	OutList1.append(Ind+"\t"+OutDict1[Ind][:-1]+"\n")
	OutList2.append(Ind+"\t"+OutDict2[Ind][:-1]+"\n")
OutFileWriting(OutFileName1, OutList1)
OutFileWriting(OutFileName2, OutList2)

#printing information on which paralogs have a certain number of individuals
#IndsPerGroupDict[Locus][ParName][GroupName]
IndsPerParDict = defaultdict(dict)#IndsPerParDict[GroupName][NumInds] = list of ParNames
GroupsPerParDict = defaultdict(list)#GroupsPerParDict[NumGroups] = list of ParNames
for Locus in IndsPerGroupDict:
	for ParName in IndsPerGroupDict[Locus]:
		GroupsPerParDict[len(IndsPerGroupDict[Locus][ParName])-1].append(ParName)
		for GroupName in IndsPerGroupDict[Locus][ParName]:
			NumInds = IndsPerGroupDict[Locus][ParName][GroupName]
			try:
				IndsPerParDict[GroupName][NumInds].append(ParName)
			except KeyError:
				IndsPerParDict[GroupName][NumInds] = [ParName]
GroupList = sorted(IndsPerParDict.keys())
GroupList.remove('all')
#first, for all individuals combined
OutFileName = OutFolder+OutFilePre+"Total_Paralogs_with_NumInds.txt"
OutList = [ ]
for NumInds in range(max(IndsPerParDict['all'].keys()), 0, -1):
	if NumInds in IndsPerParDict['all']:
		OutList.append(str(NumInds)+"\t"+",".join(IndsPerParDict["all"][NumInds])+"\n")
OutFileWriting(OutFileName, OutList)
#now doing the same for each group individually
for GroupName in GroupList:
	OutFileName = OutFolder+OutFilePre+GroupName+"_Paralogs_with_NumInds.txt"
	OutList = [ ]
	for NumInds in range(max(IndsPerParDict[GroupName].keys()), 0, -1):
		if NumInds in IndsPerParDict[GroupName]:
			OutList.append(str(NumInds)+"\t"+",".join(IndsPerParDict[GroupName][NumInds])+"\n")
	OutFileWriting(OutFileName, OutList)
#then for the number of groups, not the number of individuals
OutFileName = OutFolder+OutFilePre+"Total_Paralogs_with_NumGroups.txt"
OutList = [ ]
for NumGroups in range(max(GroupsPerParDict.keys()), 0, -1):
	if NumGroups in GroupsPerParDict:
		OutList.append(str(NumGroups)+"\t"+",".join(GroupsPerParDict[NumGroups])+"\n")
OutFileWriting(OutFileName, OutList)


#printing the IndsPerGroupDict (number of sequences of each paralog from each group)
OutFileName = OutFolder+OutFilePre+"Inds_per_Paralog.txt"
OutDict = { }
for Group in GroupList:
	OutDict[Group] = Group
KeyLine = "Group"
NumGroupsLine = "Number_of_Groups_with_Sequences"
NumSeqsLine = "Total_Sequences_in_All_Groups"
for Locus in sorted(IndsPerGroupDict.keys()):
	for ParName in sorted(IndsPerGroupDict[Locus].keys()):
		KeyLine += "\t"+ParName
		NumSeqsLine += "\t"+str(IndsPerGroupDict[Locus][ParName]['all'])
		NumGroupsLine += "\t"+str(len(IndsPerGroupDict[Locus][ParName].keys())-1)
		for Group in GroupList:
			try:
				OutDict[Group] += "\t"+str(IndsPerGroupDict[Locus][ParName][Group])
			except KeyError:
				OutDict[Group] += "\t"
OutList = [KeyLine+"\n"]
for Group in GroupList:
	OutList.append(OutDict[Group]+"\n")
OutList.append(NumSeqsLine+"\n")
OutList.append(NumGroupsLine+"\n")
OutFileWriting(OutFileName, OutList)

#writing a script to re-align the sequences to prepare them to be made into species trees!
Script1FileName = OutFolder+OutFilePre+"analysis_script1.sh"
Script2FileName = OutFolder+OutFilePre+"analysis_script2.sh"
OutScript1 = ['#! /bin/bash\n\n']
OutScript2 = [ ]
for Locus in NewParInfoDict:
	for DupName in NewParInfoDict[Locus]:
		ParName = NewParInfoDict[Locus][DupName]['Name']
		#remove existing alignment
		Line = "rm "+OutFolder+OutFilePre+ParName+OutFilePost+"_al.fa\n"
		OutScript1.append(Line)
		#align
		Line = "mafft --localpair --maxiterate 1000 --quiet "+OutFolder+OutFilePre+ParName+OutFilePost+".fa > "+OutFolder+OutFilePre+ParName+OutFilePost+"_al.fa\n"
		OutScript2.append(Line)
Line = "cat "+Script2FileName+" | parallel --joblog "+OutFolder+OutFilePre+"parallel_log.log\n"
OutScript1.append(Line)
OutFileWriting(Script1FileName, OutScript1)
OutFileWriting(Script2FileName, OutScript2)
