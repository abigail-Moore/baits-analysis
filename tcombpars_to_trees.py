#! /usr/bin/env python

#tcombpars_to_trees.py version 1.0 18 Feb. 2016 Abby Moore
#This script takes the alignments produced by tparcomb_final.py or tambig_seq_renaming.py and
#makes alignments of the subset of the sequences that are present in either the species tree
#or a list of taxa.  Then it writes a shell script to align the sequence files, build trees
#with them using RAxML, and run Notung on those trees.

import sys
from collections import defaultdict
import dendropy
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

'''
tcombpars_to_trees.py LocusListFN SpeciesTreeFN AlFolder1 AlFilePre1 AlFilePost1 AlFolder2 AlFilePre2 AlFilePost2 OutFolder OutFilePre OutFilePost ScriptPath   
'''

Usage = '''
tcombpars_to_trees.py version 1.0
This script takes the set of sequences we want to analyze further from gene
family alignments.  It writes a script to align, build trees from, and analyze 
those gene families.
tcombpars_to_trees.py
[file containing the list of loci]
[the species tree]
[first alignment folder--This one will be tried first and the second one will
not be tried if an alignment is found here!]
[prefix for first set of alignments, or "none", if none]
[suffix for first set of alignments, or "none", if none]
[second alignment folder, or "same", if no second folder]
[prefix for second set of alignments, or "none", if none]
[suffix for second set of alignments, or "none", if none]
[output folder, or "same" if the same as the alignment folder]
[prefix for the output files, or "none", if none]
[suffix for the output files, or "none", if none]
[path to the scripts, or "none", if they are on the default path]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 13:
	sys.exit("ERROR!  This script requires 12 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
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
AlFolder2 = sys.argv[6]
if AlFolder2 == "same":
	AlFolder2 = AlFolder1
elif AlFolder2[-1] != "/":
	AlFolder2 += "/"
AlFilePre2 = sys.argv[7]
if AlFilePre2 == "none":
	AlFilePre2 = ""
AlFilePost2 = sys.argv[8]
if AlFilePost2 == "none":
	AlFilePost2 = ""
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

##############################################################################################################################

#reading the list of loci
LocusList = CaptureColumn(LocusListFN, 0)
print("The list of %d loci was read from file %s.\n" % (len(LocusList), LocusListFN))
sys.stderr.write("The list of %d loci was read from file %s.\n" % (len(LocusList), LocusListFN))

#reading the list of species from the species tree
#IndList = IndListfromTree(SpeciesTreeFN)

#reading the list of species and finding the order for the outgroups from the species tree
(OGIndDict, NumInds) = OGDictfromTree(SpeciesTreeFN)
print("The names of %d individuals were read from the tree %s.\n" % (NumInds, SpeciesTreeFN))
sys.stderr.write("The names of %d individuals were read from the tree %s.\n" % (NumInds, SpeciesTreeFN))

#reading the old sequence files and writing new ones that just have the individuals we want to analyze further
OutGroupDict = { }
LocusListRemove = [ ]
NumSeqsDict = defaultdict(list)
for Locus in LocusList:
	FileExists = True
	InFileName1 = AlFolder1+AlFilePre1+Locus+AlFilePost1
	InFileName2 = AlFolder2+AlFilePre2+Locus+AlFilePost2
	try:
		SeqDictIn = SeqFileReading(InFileName1, 'fasta')
		if Vociferous == True: print("%d sequences were read for locus %s from file %s.\n" % (len(SeqDictIn), Locus, InFileName1))
	except IOError:
		try:
			SeqDictIn = SeqFileReading(InFileName2, 'fasta')
			if Vociferous == True: print("%d sequences were read for locus %s from file %s.\n" % (len(SeqDictIn), Locus, InFileName2))
		except IOError:
			LocusListRemove.append(Locus)
			print("Locus %s removed from the list." % (Locus))
			FileExists = False
	if FileExists:
		SeqDictOut = { }
		OGNum = NumInds+20
		for SeqName in SeqDictIn:
			try:
				OGRank = OGIndDict[GetIndName(SeqName)]
				SeqDictOut[SeqName] = SeqDictIn[SeqName]
				if OGRank < OGNum:
					OGNum = OGRank
					LocusOG = SeqName
			except KeyError:
				"We do not want this sequence."
		if Vociferous == True:
			print("%d of these sequences belong to the individuals of interest.\n" % (len(SeqDictOut)))
			print("The outgroup sequence for Locus %s is %s, with a rank of %d.\n" % (Locus, LocusOG, OGNum))
		NumSeqsDict[len(SeqDictOut)].append(Locus)
		OutGroupDict[Locus] = LocusOG
		OutFileName = OutFolder+OutFilePre+Locus+OutFilePost+".fa"
		SeqFileWriting(OutFileName, SeqDictOut, 'fasta')
print("New sequence files for %d loci were written, with names such as %s.\n" % (len(LocusList), OutFileName))
sys.stderr.write("New sequence files for %d loci were written, with names such as %s.\n" % (len(LocusList), OutFileName))

#re-ordering the LocusList so the files are in order from largest to smallest
LocusList = [ ]
for NumSeqs in sorted(NumSeqsDict.keys(), reverse=True):
	if NumSeqs != 0:
		LocusList += NumSeqsDict[NumSeqs]

#writing the script to analyze the data further
Script1FileName = OutFolder+OutFilePre+"analysis_script1.sh"
Script2FileName = OutFolder+OutFilePre+"analysis_script2.sh"
Script3FileName = OutFolder+OutFilePre+"analysis_script3.sh"
OutList1 = ["#! /bin/bash\n\n"]
OutList2 = []
OutList3 = []
for Locus in LocusList:
	#Align, change to phylip
	Line = "rm "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa "+OutFolder+"RAxML_*."+OutFilePre+Locus+OutFilePost+"\n"
	OutList1.append(Line)
	Line = "mafft --localpair --maxiterate 1000 --quiet "+OutFolder+OutFilePre+Locus+OutFilePost+".fa > "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa && "
	Line += ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.fa\n"
	OutList2.append(Line)
	#raxml, Notung
	Line = "raxmlHPC -f a -s "+OutFolder+OutFilePre+Locus+OutFilePost+"_al.phy -n "+OutFilePre+Locus+OutFilePost+" -m GTRCAT -p 1234 -N 100 -x 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+" && "
	#version for laptop
	#Line += "java -jar ~/bin/Notung-2.8.1.6-beta.jar -g "+OutFolder+"RAxML_bipartitions."+OutFilePre+Locus+OutFilePost+" -s "+SpeciesTreeFN+" --speciestag prefix --rearrange --threshold 90 --homologtabletabs --nolosses --silent  --usegenedir\n"
	#version for oscar
	Line += "java -jar "+ScriptPath+"Notung-2.8.1.6-beta.jar -g "+OutFolder+"RAxML_bipartitions."+OutFilePre+Locus+OutFilePost+" -s "+SpeciesTreeFN+" --speciestag prefix --rearrange --threshold 90 --homologtabletabs --nolosses --silent  --usegenedir\n"
	OutList3.append(Line)
Line = "cat "+Script2FileName+" | parallel --joblog "+OutFolder+OutFilePre+"parallel_log2.log\n"
Line += "cat "+Script3FileName+" | parallel --joblog "+OutFolder+OutFilePre+"parallel_log3.log\n"
OutList1.append(Line)
OutFileWriting(Script1FileName, OutList1)
OutFileWriting(Script2FileName, OutList2)
OutFileWriting(Script3FileName, OutList3)
