#! /usr/bin/env python

#tclade_finder.py version 1.0 8 May 2015 Abby Moore
#This script is for finding clades of sequences that come from the baits using Dendropy
#to read the trees.  Then it determines whether the clade is supported as belonging to
#a particular paralog of the gene and makes a sequence alignment of that paralog (including
#both the sequences of interest and the previous sequences we had for the paralog).

import dendropy #This is obviously going to be important
from collections import defaultdict
import sys
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

#examples:
'''
tclade_finder.py ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/seqfiles_head.txt /home/abby/transcriptomes/baits_Bv/pardict3.txt ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/ of3t_ same pruned_RAxML_bipartitions.of3_ ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/ppc_contigs/ of4_ ~/transcriptomes/C4alignmentsPA/pepc/ none none
tclade_finder.py ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/seqfiles_head2.txt /home/abby/transcriptomes/baits_Bv/pardict3.txt ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/ of3_ same reroot_RAxML_bipartitions.of3_ ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/ppc_contigs/ of4_ ~/transcriptomes/C4alignmentsPA/pepc/ none none Lewisia 
tclade_finder.py SeqFileList PFileName SeqFolder SeqFilePre TreeFolder TreeFilePre OutFolder OutFilePre AlFolder AlFilePre AlFilePost InGroupList
'''

Usage = '''
tclade_finder.py is a script to find monophyletic clades that consist entirely 
of sequences of interest and determine the sequences in those clades.
IMPORTANT!! If one of the original (backbone) sequences is nested within a clade
formed by the sequences being analyzed, it will break those clades up into very
many smaller clades.  In this case, this sequence should be added to the 
alignment of sequences of interest and it will be analyzed with them.
tclade_finder.py 
[list of sequence files (in the first column; there can be multiple columns, but
it will only read the first one)]
[file showing which bait sequences belong to which paralogs]
[folder where sequence files are found] 
[prefix for sequence files, or "none", if none]
[folder where tree files are found or "same" if the same as the sequence folder]
[prefix for tree files, or "none", if none]
[output folder, or "same", if it is the same as the sequence files are found]
[prefix for output files, or "none", if none]
[for phylogenetic analysis: folder containing the template alignments]
[template file prefix or "none"]
[template file ending or "none"]
[list of genera in ingroup, so these sequences can be skipped]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 15:
	sys.exit("ERROR! This script requires at least 11 additional arguments, and you supplied %d.  %s" % (len(sys.argv)-1, Usage))
SeqFileList = sys.argv[1]
PFileName = sys.argv[2]
SeqFolder = sys.argv[3]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
SeqFilePre = sys.argv[4]
if SeqFilePre == "none":
	SeqFilePre = ""
TreeFolder = sys.argv[5]
if TreeFolder == "same":
	TreeFolder = SeqFolder
elif TreeFolder[-1] != "/":
	TreeFolder += "/"
TreeFilePre = sys.argv[6]
OutFolder = sys.argv[7]
if OutFolder == "same":
	OutFolder = SeqFolder
elif OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[8]
if OutFilePre == "none":
	OutFilePre = ""
AlFolder = sys.argv[9]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[10]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[11]
if AlFilePost == "none":
	AlFilePost = ""
if len(sys.argv) > 12:
	InGroupList = sys.argv[12:]
else:
	InGroupList = [ ]

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
	#This is ContigFileList

#DDictFromFile makes a three-level dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
#modified from tbaits_intron_removal.py
def DDictFromFile(FileName):
	TempDict = defaultdict(dict)
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]][Line[1]] = Line[2]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1], Line[2]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1], Line[2]))
	return TempDict

#This makes
#modified from LevelDictList from tcontigs_to_paralogs.py
def PrunedDictList(InFileName, NumHeaderLines, Key2Column, ListColumn, WantedList):
	TempDict = defaultdict(dict)
	InFile = open(InFileName, 'rU')
	LineNum = 0
	#making the dictionary from the file
	for Line in InFile:
		if LineNum >= NumHeaderLines:
			Line = Line.strip('\r').strip('\n').split('\t')
			if Line[Key2Column] in WantedList:
				try:
					TempDict[Line[0]][Line[Key2Column]].append(Line[ListColumn])
				except KeyError:
					TempDict[Line[0]][Line[Key2Column]] = [Line[ListColumn]]
		LineNum += 1
	InFile.close()
	#condensing the dictionary (by removing duplicate values from the lists), if necessary
	for Key0 in TempDict:
		for Key1 in TempDict[Key0]:
			ListTemp = TempDict[Key0][Key1]
			ListTemp = list(set(ListTemp))
			TempDict[Key0][Key1] = ListTemp
	#print("%d lines were read from the file %s, making a dictionary of the form:\n\
#Key1: %s, Key2: %s, Value(list): %s.\n" % (LineNum, InFileName, Key0, Key1, ", ".join(ListTemp)))
	#sys.stderr.write("%d lines were read from the file %s, making a dictionary of the form:\n Key1: %s, Key2: %s, Value(list): %s.\n" % (LineNum, InFileName, Key0, Key1, ", ".join(ListTemp)))
	return TempDict
	#This is InGroupSeqDict

#LocusSeqGetter reads a series of sequence files and makes a dictionary of the sequences
#that have been classified according to locus.  Used first in tbaits_intron_removal.py
def LocusSeqGetter(FileList,Folder,FilePre,FilePost,SeqFormat):
	TempDict = defaultdict(dict)
	for SeqFile in FileList:
		Locus = SeqFile.replace(FilePre,"").replace(FilePost,"")
		FileName = Folder + SeqFile
		InFile = open(FileName, 'rU')
		for record in SeqIO.parse(InFile, SeqFormat):
			TempDict[Locus][record.id] = str(record.seq)
		#print("%d sequences for locus %s were read from the file %s.\n" % (len(TempDict[Locus].keys()), Locus, FileName))
	print("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	sys.stderr.write("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	return TempDict
	#This is ContigDict

#GroupFinder goes through a tree and finds monophyletic groups of the sequences of interest.
def GroupFinder(SeqDict, Folder, FilePre, TreeFormat, OtherSeqsDict):
	GroupDict = defaultdict(dict)#GroupDict[Locus][SeqName] = GroupNum--This is not exported
	TempDict = defaultdict(dict)#TempDict[Locus][GroupNum][various things] = various things--This is exported
	for Locus in SeqDict:
		SeqsWanted = SeqDict[Locus].keys()
		#making a list of the sequences we need to understand plus any other sequences from that genus that are already in the tree
		InGroupSeqList = [ ]
		#and a list of just the additional sequences, so they can be removed at the end
		AddlSeqList = [ ]
		try:
			for GenusName in OtherSeqsDict[Locus]:
				InGroupSeqList += OtherSeqsDict[Locus][GenusName]
				AddlSeqList += OtherSeqsDict[Locus][GenusName]
		except KeyError:
			"do nothing"
		InGroupSeqList += SeqDict[Locus].keys()
		TreeFileName = Folder+FilePre+Locus
		#reading the tree for this locus
		tree1 = dendropy.Tree(stream=open(TreeFileName), schema=TreeFormat, preserve_underscores=True, as_rooted=False)
		#print(tree1.as_ascii_plot(show_internal_node_labels=True))
		#print Locus
		CurrentGroup = 1
		TempDict[Locus] = defaultdict(dict)
		GoodGroups = 0
		BadGroups = 0
		#finding out where all of our sequences of interest are in the tree
		for SeqName in SeqsWanted:
			#only if that sequence hasn't already been classified
			if GroupDict[Locus].get(SeqName,'none') == 'none':
				CladeList = [SeqName]
				GroupQual = 'bad'
				#print SeqName
				node = tree1.find_node_with_taxon_label(SeqName)
				#then finding its sister taxa
				for node3 in node.sister_nodes():
					SisTaxList = [str(node2.taxon) for node2 in node3.leaf_nodes()]
				#if all of the members of the clade are in our group of interest, step back through the tree to find the whole clade
				if (all(SisTax in InGroupSeqList for SisTax in SisTaxList) == True):
					CladeList += SisTaxList
					#now, we need to go back into the tree a node at a time
					for node_c in (node.parent_node).ancestor_iter():
						nodep = node_c.parent_node
						SisTaxList = [str(node2.taxon) for node2 in nodep.leaf_nodes()]
						#continue further back if all of the descendants of that node are in the clade of interest
						if (all(SisTax in InGroupSeqList for SisTax in SisTaxList) == True):
							CladeList += SisTaxList
						#but stop if they aren't and figure out which ones are outgroups
						else:
							CladeList = list(set(CladeList))
							#remove extraneous sequences, because we don't need to classify them
							if AddlSeqList != [ ]:
								for SeqName in AddlSeqList:
									if SeqName in CladeList:
										CladeList.remove(SeqName)
							TempDict[Locus][CurrentGroup]['Mems'] = CladeList
							TempDict[Locus][CurrentGroup]['BS-group'] = node_c.label 
							TempDict[Locus][CurrentGroup]['BS-parent'] = [nodep.label]
							ParentNode = nodep
							#print ("Group %d: support %s\n" % (CurrentGroup, nodep.label))
							break
				#otherwise, this sequence forms a clade by itself
				else:
					TempDict[Locus][CurrentGroup]['SisterGroup'] = SisTaxList
					TempDict[Locus][CurrentGroup]['Mems'] = [SeqName]
					TempDict[Locus][CurrentGroup]['BS-group'] = 100 #not really, but there is only one member
					TempDict[Locus][CurrentGroup]['BS-parent'] = [(node.parent_node).label]
					ParentNode = node.parent_node
				#then step through the tree to find a well-supported group that includes this clade and its sister group, to determine which paralog we have
				if int(TempDict[Locus][CurrentGroup]['BS-parent'][0]) < 75:
					for node_d in (ParentNode.parent_node).ancestor_iter():
						TempDict[Locus][CurrentGroup]['BS-parent'].append(node_d.label)
						SisTaxList = [str(node2.taxon) for node2 in node_d.leaf_nodes()]#I hope that by moving this up here, this will work, even when it goes to the bottom of the tree without finding a well-supported node.
						#first, escape if we are at the bottom of the tree
						if (node_d.label == None):
							break
						#continue if the node is not well-supported
						elif int(node_d.label) < 75:
							"do nothing"
							#print ("Group %d, support %s, continue\n" % (CurrentGroup, node_d.label))
						#finish if the node is well-supported
						else:
						#if int(node_d.label) > 75:
							GroupQual = 'good'
							#print ("Group %d, support %s, break\n" % (CurrentGroup, node_d.label))
							#SisTaxList = [str(node2.taxon) for node2 in node_d.leaf_nodes()] #I don't think I need this line anymore, because I have it above as well.
							break
				else:
					GroupQual = 'good'
				if GroupQual == 'good':
					GoodGroups += 1
				elif GroupQual == 'bad':
					BadGroups += 1
				OGList = [ ]
				for Name in SisTaxList:
					if (Name in CladeList) == False:
						OGList.append(Name)
				TempDict[Locus][CurrentGroup]['SisterGroup'] = OGList
				for ContigName in CladeList:
					GroupDict[Locus][ContigName] = CurrentGroup
				CurrentGroup += 1
		print("A total of %d monophyletic groups of sequences was found for locus %s.\n" % (CurrentGroup-1, Locus))
		sys.stderr.write("A total of %d monophyletic groups of sequences was found for locus %s.\n" % (CurrentGroup-1, Locus))
		print("%d of these formed well-supported groups with backbone sequences and %d did not.\n" % (GoodGroups, BadGroups))
		sys.stderr.write("%d of these formed well-supported groups with backbone sequences and %d did not.\n" % (GoodGroups, BadGroups))
		if CurrentGroup-1 > 5:
			print("There are quite a large number of groups (%d) for locus %s.  Perhaps one of the backbone sequences is actually part of the clade of interest and it should be reclassified?\n" % (CurrentGroup-1, Locus))
			sys.stderr.write("There are quite a large number of groups (%d) for locus %s.  Perhaps one of the backbone sequences is actually part of the clade of interest and it should be reclassified?\n" % (CurrentGroup-1, Locus))
	return TempDict
	#This is SisterGroupDict

def GroupParalogFinder(PDict, SGDict, SeqDict):
	#SGDict[Locus][GroupNum]['SisterGroup'] = [list of sister group members]
	#SGDict[Locus][GroupNum]['Mems'] = [list of group members]
	#SGDict[Locus][GroupNum]['BS-parent'] = [list of BS values for subtending nodes (increasingly deeper in tree)]
	#SGDict[Locus][GroupNum]['BS-group'] = BS value for clade
	#PDict[Locus][SeqName] = Paralog
	for Locus in SGDict:
		ContigList = SeqDict[Locus].keys()
		WithParalogs = 0
		AmbigSeqs = 0
		for GroupNum in SGDict[Locus]:
			PList = [ ]
			for SeqName in SGDict[Locus][GroupNum]['SisterGroup']:
				try:
					PList.append(PDict[Locus][SeqName])
				except KeyError:
					if ((SeqName in ContigList) == False) and SeqName[0:2] != "Bv":
						print("WARNING! Sequence %s was not in the list of template sequences that were classified into paralogs and was also not a contig!\n" % (SeqName))
						sys.stderr.write("WARNING! Sequence %s was not in the list of template sequences that were classified into paralogs and was also not a contig!\n" % (SeqName))
			PList = list(set(PList))
			if len(PList) == 1:
				GParalog = PList[0]
				WithParalogs += 1
			else:
				GParalog = Locus+'-Ambig'+str(GroupNum)
				AmbigSeqs += 1
			if (SGDict[Locus][GroupNum]['BS-parent'][-1] == None) or (int(SGDict[Locus][GroupNum]['BS-parent'][-1]) < 75):
				GParalog = Locus+'-Ambig'+str(GroupNum)
				AmbigSeqs += 1
			SGDict[Locus][GroupNum]['Paralog'] = GParalog
			#print("%s: %d, %s, %s, %s\n" % (Locus, GroupNum, " ".join(SGDict[Locus][GroupNum]['BS-parent']), " ".join(PList), GParalog)) 
		print("Of the clades for locus %s, %d could be identified to paralog, while %d could not be.\n" % (Locus, WithParalogs, AmbigSeqs))
		sys.stderr.write("Of the clades for locus %s, %d could be identified to paralog, while %d could not be.\n" % (Locus, WithParalogs, AmbigSeqs))
	return SGDict
	#This is SisterGroupDict2

#SeqstoParalogs combines the results of LocusSeqGetter, GroupFinder, and GroupParalogFinder
#to make a dictionary where the sequences are classified first by locus and then by
#paralog.
def SeqstoParalogs(CDict, SGDict):
	#CDict[Locus][Contig] = Seq
	#SGDict[Locus][GroupNum][Various] = Various
	TempDict = defaultdict(dict)#TempDict[Locus][Paralog][ContigName] = Sequence
	for Locus in SGDict:
		TempDict[Locus] = defaultdict(dict)
		for GroupNum in SGDict[Locus]:
			Paralog = SGDict[Locus][GroupNum]['Paralog']
			for Contig in SGDict[Locus][GroupNum]['Mems']:
				ContigName = Contig+"-"+str(GroupNum)
				TempDict[Locus][Paralog][ContigName] = CDict[Locus][Contig]
			#print("%s: %s: %d: %d" % (Locus, Paralog, len(TempDict[Locus][Paralog]), GroupNum))
	return TempDict
	#This is OutDict

#LocusParalogSeqWriter writes sequences in the desired format to files from a
#dictionary where the top-level key is the locus, the second-level key is the group number,
#the third-level key is the sequence name, and the value is the sequence.
#The top-level key can simply be "".
#modified from LocusGroupSeqWriter from tbaits_intron_removal.py
def LocusParalogSeqWriter(SeqDict, Folder, Prefix, Suffix, SeqFormat):
	#SeqDict[Locus][Paralog][ContigName] = Sequence
	TempDict = defaultdict(dict)#TempDict[Locus][Paralog] = file name
	for Locus in SeqDict:
		NumFiles = 0
		for Paralog in SeqDict[Locus]:
			OutFileName = Folder+Prefix+Paralog+Suffix
			TempDict[Locus][Paralog] = Prefix+Paralog
			OutFile = open(OutFileName, 'w')
			for ContigName in SeqDict[Locus][Paralog]:
				Record1 = SeqRecord(seq=Seq(SeqDict[Locus][Paralog][ContigName], IUPAC), id = ContigName, description = "")
				SeqIO.write(Record1, OutFile, SeqFormat)
			OutFile.close()
			NumFiles += 1
		print("%d sequence files were written for the locus %s, with names such as %s.\n" % (NumFiles, Locus, OutFileName))
		sys.stderr.write("%d sequence files were written for the locus %s, with names such as %s.\n" % (NumFiles, Locus, OutFileName))
	return TempDict
	#This is SeqFileDict
	
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

#This is based on MRScriptWriter from tbaits_introns_removal.py
#It writes a script that adds the sequences to the file of pre-existing
#sequences for that locus and then aligns them.
def AMScriptWriter(SeqFileDict, Folder, Prefix, AFolder, APre, APost):
	OutList = ["#! /bin/bash\n"]
	for Locus in SeqFileDict:
		for Paralog in SeqFileDict[Locus]:
			NamePart = SeqFileDict[Locus][Paralog]
			Line = "rm "+Folder+NamePart+"_pars_al.fa\n"
			Line += "rm "+Folder+"RAxML*"+NamePart+"\n"
			if Paralog[:len(Locus+"-Ambig")] != Locus+"-Ambig":
				Line += "mafft --addfragments "+Folder+NamePart+".fa --quiet --thread -1 "+AFolder+APre+Paralog+APost+".fa > "+Folder+NamePart+"_pars_al.fa\n"
			else:
				Line += "mafft --addfragments "+Folder+NamePart+".fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+NamePart+"_pars_al.fa\n"
			OutList.append(Line)
	OutFileName = Folder+Prefix+"Analysis_Script.sh"
	OutFileWriting(OutFileName, OutList)
	print("The shell script for analyzing the sequences further was written to %s.\n" % (OutFileName))
	sys.stderr.write("The shell script for analyzing the sequences further was written to %s.\n" % (OutFileName))
	
#read the names of the sequence files from blast file list read in tbaits_intron_removal.py
ContigFileList = CaptureColumn(SeqFileList, 0)

#make a dictionary of which backbone sequences belong to which paralogs
if PFileName != "none":
	ParalogDict = DDictFromFile(PFileName)
else:
	print("No file was given for classifying the backbone sequences into paralogs.  It is assumed that their classification is unknown.\n")
	sys.stderr.write("No file was given for classifying the backbone sequences into paralogs.  It is assumed that their classification is unknown.\n")

#make a dictionary of the other sequences from the genus/genera of interest that are already in the tree
InGroupSeqDict = PrunedDictList(PFileName, 0, 3, 1, InGroupList)

#read the sequence files
ContigDict = LocusSeqGetter(ContigFileList, SeqFolder, SeqFilePre, '.fa', 'fasta')

#finding groups in the tree files
#SisterGroupDict = GroupFinder(ContigDict, TreeFolder, TreeFilePre, 'nexus', InGroupSeqDict)
SisterGroupDict = GroupFinder(ContigDict, TreeFolder, TreeFilePre, 'newick', InGroupSeqDict)

#figuring out which groups really have good classifications
SisterGroupDict2 = GroupParalogFinder(ParalogDict, SisterGroupDict, ContigDict)

#write the output
#sequence files for each group/locus:

OutDict = SeqstoParalogs(ContigDict, SisterGroupDict2)
FileDict = LocusParalogSeqWriter(OutDict, OutFolder, OutFilePre, '.fa', 'fasta')

#script for analyzing the sequences
AMScriptWriter(FileDict, OutFolder, OutFilePre, AlFolder, AlFilePre, AlFilePost)

#an information file:
#Locus\tGroupNum\tBS_Group\tBS_Clade\tNumNodes (that I had to go down to get that BS value)\tParalog
#another information file:
#Locus\tContig\tGroupNum\tParalog
OutList = ['Locus\tGroupNum\tBS_Group\tBS_Clade\tNumSisterNodes\tParalog\n' ]
OutGroupList = ['Locus\tContig\tGroupNum\tParalog\n']
for Locus in SisterGroupDict2:
	for GroupNum in SisterGroupDict2[Locus]:
		Line = Locus+"\t"+str(GroupNum)+"\t"+str(SisterGroupDict2[Locus][GroupNum]['BS-group'])+"\t"+str(SisterGroupDict2[Locus][GroupNum]['BS-parent'][-1])+"\t"
		Line += str(len(SisterGroupDict2[Locus][GroupNum]['BS-parent']))+"\t"+SisterGroupDict2[Locus][GroupNum]['Paralog']+"\n"
		OutList.append(Line)
		for Contig in SisterGroupDict2[Locus][GroupNum]['Mems']:
			LineG = Locus+"\t"+Contig+"-"+str(GroupNum)+"\t"+str(GroupNum)+"\t"+SisterGroupDict2[Locus][GroupNum]['Paralog']+"\n"
			OutGroupList.append(LineG)
OutFileName = OutFolder+OutFilePre+"Paralog_Info.txt"
OutFileGName = OutFolder+OutFilePre+"Contig_Groups.txt"
OutFileWriting(OutFileName, OutList)
OutFileWriting(OutFileGName, OutGroupList)
