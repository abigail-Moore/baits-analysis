#! /usr/bin/env python

#tnotung_homolog_parsing.py version 1.0 2 Feb. 2016
#This script parses the trees produced by Notung and determines which sequences are orthologous and which are paralogous.
#They need to be rooted and in Notung format, e.g., RAxML_bipartitions.ppt_nsbest.rearrange.0.ntg

import sys
import re
from collections import defaultdict, OrderedDict
import dendropy
from biolite import workflows
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
[output folder]
[prefix for output files]
[suffix for output files--sequence files only]

For all but the first folder name, substituting "same" causes them to be the
same as the first folder.
For all prefixes and suffixes, substituting "none" causes there to be no 
prefix/suffix.
'''

'''
tnotung_homolog_parsing.py InFolder InFilePre InFilePost AlFolder AlFilePre AlFilePost LocusListFileName OutFolder OutFilePre OutFilePost
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 11:
	sys.exit("ERROR!  This script requires 12 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
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
OutFolder = sys.argv[8]
if OutFolder ==  "same":
	OutFolder = InFolder
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[9]
if OutFilePre == "none":
	OutFilePre = ""
OutFilePost = sys.argv[10]
if OutFilePost == "none":
	OutFilePost = ""

ShowTrees = False
Vociferous = False

#####################################################################################

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
	TreeFile = open(FileName, 'rU')
	LineNum = 1
	for Line in TreeFile:
		#The tree is in the first line.
		if LineNum == 1:
			#Various regular expressions to get rid of stuff we don t need and to make the correct node labels.
			r1 = re.compile(r"[rn]\d+:(\d+\.\d+)\[&&NHX:S=\w+:Nset=<[\w@]+>:H=\w:D=(\w)[:B=0-9\.\]]+")
			r3= re.compile(r"[rn]\d+\[&&NHX:S=\w+:Nset=<[\w@]+>:H=\w:D=(\w)\]")
			r5 = re.compile(r"\[&&NHX:S=\w+:Nset=<[\w@]+>:H=\w+\]")
			TreeString = r5.sub(r"",r3.sub(r"\1",r1.sub(r"\2:\1",Line.rstrip())))
		LineNum += 1
	TreeFile.close()
	#turning the string into a tree
	TempTree = dendropy.Tree.get_from_string(TreeString, schema='newick', preserve_underscores=True)
	if ShowTrees == True: print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	return TempTree
	#This is LocusTree_orig

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
	OutSeqDict = { }#OutSeqDict[SeqName] = Seq
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
					(ConSeq, NumAmbig, AmbigSitesList, NumSeqs, Overlap) = ConSeqMaker(InSeqDict, DaughterDict[IndName])
					if NumAmbig < 5:
						CombinableInds.append(IndName)
						if Vociferous == True: print("The sequences for individual %s were combinable!\n" % (IndName))
						ParNameList =sorted(list(set([GetParName(SeqName) for SeqName in DaughterDict[IndName]])))
						#Getting the NewParName, and trying to avoid Ambig, if possible
						NewParName = ""
						for ParName in ParNameList:
							if ParName[:5] != "Ambig":
								NewParName = ParName
						#but, if all of the potential paralogs are Ambig, then choosing the first one
						if NewParName == "":
							NewParName = ParNameList[0]
						NewSeqName = ".".join([IndName, ParNameList[0], str(NumSeqs)+"seqs", str(NumAmbig)+"ambig", "len"+str(len(ConSeq))])
						OutSeqDict[NewSeqName] = ConSeq
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
					else:
						if Vociferous == True: print("The sequences for individual %s were not combinable.  :(\n" % (IndName))
			#going through the leaf nodes again and removing/renaming them, if necessary
			for Daughter_Node in Node.leaf_nodes():
				OldName = Daughter_Node.taxon.label
				if OldName in RenameDict.keys():
					Daughter_Node.taxon.label = RenameDict[OldName]
				elif OldName in RemoveList:
					TempTree.prune_taxa([Daughter_Node.taxon], update_bipartitions=False)
			#relabeling the node
			if sorted(CombinableInds) == sorted(DupInds):
				if Vociferous == True: print("Duplication %s successfully resolved!\n" % (Node.label))
				Node.label = "N"
				NumDupsNew = NumDupsNew -1
	if ShowTrees == True: print(TempTree.as_ascii_plot(show_internal_node_labels=True))
	if Vociferous == True: print("Prior to combining separated parts of sequences of the same paralog, there were %d duplications.  Now there are %d.\n" % (NumDupsOld, NumDupsNew))
	return (TempTree, OutSeqDict)

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
						NewParName = "?"+"_".join(sorted(DictTemp[DupNum][ParNum]['ParNames']))
						DictTemp[DupNum][ParNum]['NewParNames'].append(NewParName)
						for SeqName in DictTemp[DupNum][ParNum]['Seqs']:
							OutDictTemp[SeqName] = NewParName
						if Vociferous == True: print("The new paralog name for half %d of duplication number %s is uncertain, but is temporarily %s.\n" % (ParNum, DupNum, NewParName))
	#Then going through all of the SeqNames in the tree to make sure they have an associated paralog in the dictionary
	for LeafNode in TempTree.leaf_nodes():
		SeqName = LeafNode.taxon.label
		#possibly the sequence was given a new name starting with "?"
		try:
			if OutDictTemp[SeqName][0] == "?":
				print("ERROR!  Sequence %s was not properly classified according to paralog!\n" % (SeqName))
				OutDictTemp[SeqName] = OutDictTemp[SeqName][1:]
		#If not, they were not duplicated, and their paralog should be the locus name.
		except KeyError:
			OutDictTemp[SeqName] = Locus
	return DictTemp, OutDictTemp, RDict
	#These are ParInfoDict[Locus], NewParDict[Locus], and RankDict.

#OrthoParFinder goes through a tree and separates it into sets of orthologous genes, while making separate groupings for the
#various paralogs
#original
def OrthoParFinder(TempTree, PNDict, RDict):
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
	if Vociferous == True: print("The new name for this paralog will be %s.\n" % (BestPar))
	DictTemp['root']['Name'] = BestPar
	return DictTemp

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
	#dumb_consensus works well as long as there are no ambiguities.  If I want to be able to notice/count them, I need something more sophisticated.
	#AlignTempInfo = AlignInfo.SummaryInfo(AlignTemp)
	#ConSeq = AlignTempInfo.dumb_consensus(ambiguous='-', consensus_alpha=IUPAC)
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
	return (ConSeq, AmbigNucs, AmbigNucList, NSeqs, Overlap)

#######################################################################################

LocusList = CaptureColumn(LocusListFileName, 0)

#The dictionaries that will be filled out:
TreeDict = { }#TreeDict[Locus] = LocusTree--nodes labeled according to duplications, with spurious duplications removed
ParInfoDict = defaultdict(dict)#ParInfoDict[Locus][DupNum][ParNum][various keys] = various lists
NewParDict = defaultdict(dict)#NewParDict[Locus][SeqName] = ParName
NewParInfoDict = defaultdict(dict)#NewParInfoDict[Locus][DupName]['SeqNames'/'ParList'/'Name'] = list of SeqNames/list of ParNames/final name of paralog
DupPosDict = defaultdict(dict)#DupPosDict[Locus][DupName]['sister'/'daughters'] = list of individuals

for Locus in LocusList:
	#reading the tree from Notung (.ntg format)
	TreeFileName = InFolder+"RAxML_bipartitions."+InFilePre+Locus+InFilePost+".rearrange.0.ntg"
	LocusTree_orig = NotungFileReading(TreeFileName)
	#reading the alignment (fasta format)
	AlFileName = AlFolder+AlFilePre+Locus+AlFilePost+".fa"
	LocusAl_orig = SeqFileReading(AlFileName, 'fasta')
	#going through each node in the tree and determining whether Notung reconstructed a duplication at that node
	(LocusTree_new) = DuplicationFinder(LocusTree_orig)
	#Go through the duplications and try to combine the sequences that can be combined
	(LocusTree_pruned, NewSeqDict) = SeqCombiner(LocusTree_new, LocusAl_orig)
	#then figuring out where all the duplications are
	DupPosDict[Locus] = DuplicationMapper(LocusTree_pruned)
	#choose the longer sequence for each duplication involving only one individual that could not be successfully resolved
	(LocusTree_nodups) = IndDupsRemoving(LocusTree_pruned, LocusAl_orig)
	#finding a name for each duplication
	(ParInfoDict[Locus], NewParDict[Locus], RankDict) = ParalogNamer(LocusTree_new, Locus)
	#go through and find more or less orthologous sequences, using the half of each duplication with the largest taxon sampling
	NewParInfoDict[Locus] = OrthoParFinder(LocusTree_new, NewParDict[Locus], RankDict)
	#write a new sequence file for the gene family
	OutFileNameGF = OutFolder+OutFilePre+Locus+OutFilePost+".fa"
	SeqDictGF = { }
	#write new sequence files for each new paralog
	for DupNode in NewParInfoDict[Locus]:
		ParName = NewParInfoDict[Locus][DupNode]['Name']
		OutFileName = OutFolder+OutFilePre+Locus+"_"+ParName+OutFilePost+".fa"
		SeqDict = { }
		for SeqName in NewParInfoDict[Locus][DupNode]['SeqNames']:
			#first assuming the sequence was an existing sequence
			try:
				SeqDict[SeqName] = LocusAl_orig[SeqName]
				SeqDictGF[SeqName] = LocusAl_orig[SeqName]
			#if not, it must be a new combined sequence
			except KeyError:
				SeqDict[SeqName] = NewSeqDict[SeqName]
				SeqDictGF[SeqName] = NewSeqDict[SeqName]
		SeqFileWriting(OutFileName, SeqDict, 'fasta')
	SeqFileWriting(OutFileNameGF, SeqDictGF, 'fasta')

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

#writing a script to re-align the sequences to prepare them to be made into species trees!
Script1FileName = OutFolder+OutFilePre+"analysis_script1.sh"
Script2FileName = OutFolder+OutFilePre+"analysis_script2.sh"
OutScript1 = ['#! /bin/bash\n\n']
OutScript2 = [ ]
for Locus in NewParInfoDict:
	for DupName in NewParInfoDict[Locus]:
		ParName = NewParInfoDict[Locus][DupNode]['Name']
		#remove existing alignment
		Line = "rm "+OutFolder+OutFilePre+Locus+"_"+ParName+OutFilePost+"_al.fa\n"
		OutScript1.append(Line)
		#align
		Line = "mafft --localpair --maxiterate 1000 --quiet "+OutFolder+OutFilePre+Locus+"_"+ParName+OutFilePost+".fa > "+OutFolder+OutFilePre+Locus+"_"+ParName+OutFilePost+"_al.fa\n"
		OutScript2.append(Line)
Line = "cat "+Script2FileName+" | parallel --joblog "+OutFolder+OutFilePre+"parallel_log.log\n"
OutScript1.append(Line)
OutFileWriting(Script1FileName, OutScript1)
OutFileWriting(Script2FileName, OutScript2)
