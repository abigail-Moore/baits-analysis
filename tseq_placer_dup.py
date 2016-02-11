#! /usr/bin/env python

#tseq_placer_dup.py version 1.0 19 June 2015 Abby Moore
#This script reads the output from RAxML's classification of short reads and determines whether or not
#each contig can be accurately classified according to paralog.
#Everything but reading the tree and placing the sequences is based on tclade_finder.py
#sequence file names should be of the form: SeqFolder+SeqFilePre + Locus + SeqFilePost + ".fa"
#tree file names should be of the form: TreeFolder+"RAxML_originalLabelledTree."+TreeFilePre+Locus
#likelihood file names should be of the form: TreeFolder+"RAxML_classificationLikelihoodWeights."+TreeFilePre+Locus
#This file outputs fasta files of the contigs classified according to paralog (or into groups of paralogs, if it is ambiguous).
#Both ambiguously-classified contigs and contigs that could be classified without ambiguity are treated in the same way.
#The shell script it outputs simply aligns the contigs back to the same backbone alignment to which they were originally aligned, not
#to the separate alignments for each paralog.  (Perhaps this should be changed?**)

#tseq_placer_dup.py LocusFileList PFileName SeqFolder SeqFilePre SeqFilePost TreeFolder TreeFilePre OutFolder OutFilePre AlFolder AlFilePre AlFilePost CutOff GrpFate[merge, separate]
#tseq_placer_dup.py ~/transcriptomes/sandbox/amk1/LocusList.txt ~/transcriptomes/general/pardict6.txt ~/transcriptomes/sandbox/amk1/ Mont_ none ~/transcriptomes/sandbox/amk1/ Mont_ ~/transcriptomes/sandbox/amk1/ trial1_ ~/transcriptomes/sandbox/amk1/ none _allgoodseqs_al 0.10 merge

from collections import defaultdict
import dendropy
import re
import sys
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

Usage = '''
tseq_placer_dup.py
[list of loci (meaning the loci after which the files are named)]
[file showing which bait sequences belong to which paralogs]
[folder where sequence files are found] 
[prefix for sequence files, or "none", if none]
[suffix for sequence files, or "none", if none]
[folder where tree files are found or "same" if the same as the sequence folder]
[prefix for tree files, or "none", if none]
[output folder, or "same", if it is the same as the sequence files are found]
[prefix for output files, or "none", if none]
[for phylogenetic analysis: folder containing the template alignments]
[template file prefix or "none"]
[template file ending or "none"]
[cutoff likelihood value for a given contig placement for it to be considered 
in the classification of that contig]
[fate of groups that contain the same taxa (but do not contain overlapping 
nodes): "merge" or "separate"]
'''

GrpFateList = ['merge', 'separate']

print("%s\n" % (" ".join(sys.argv)))
sys.stderr.write("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 15:
	sys.exit("ERROR! This script requires 14 additional arguments, and you supplied %d.  %s" % (len(sys.argv)-1, Usage))
LocusFileList = sys.argv[1]
PFileName = sys.argv[2]
SeqFolder = sys.argv[3]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
SeqFilePre = sys.argv[4]
if SeqFilePre == "none":
	SeqFilePre = ""
SeqFilePost = sys.argv[5]
if SeqFilePost == "none":
	SeqFilePost = ""
TreeFolder = sys.argv[6]
if TreeFolder == "same":
	TreeFolder = SeqFolder
elif TreeFolder[-1] != "/":
	TreeFolder += "/"
TreeFilePre = sys.argv[7]
if TreeFilePre == "none":
	TreeFilePre = ""
OutFolder = sys.argv[8]
if OutFolder == "same":
	OutFolder = SeqFolder
elif OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[9]
if OutFilePre == "none":
	OutFilePre = ""
AlFolder = sys.argv[10]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[11]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[12]
if AlFilePost == "none":
	AlFilePost = ""
CutOff = float(sys.argv[13])
GrpFate = sys.argv[14]
if (GrpFate in GrpFateList) == False:
	sys.exit("ERROR!  The fate of the groups can be %s, but you wrote %s.\n %s" % (", ".join(GrpFateList), GrpFate, Usage))


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
	#sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is LocusList

#DDictFromFile makes a three-level dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
#modified from tbaits_intron_removal.py
def DDictFromFile(FileName):
	TempDict = defaultdict(dict)
	LineNum = 0
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]][Line[1]] = Line[2]
		LineNum += 1
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s: %s\n" % (LineNum, FileName, Line[0], Line[1], Line[2]))
	#sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s: %s\n" % (LineNum, FileName, Line[0], Line[1], Line[2]))
	return TempDict

#LocusSeqGetter reads a series of sequence files and makes a dictionary of the sequences
#that have been classified according to locus.  Modified from tbaits_intron_removal.py
def LocusSeqGetter(LList,Folder,FilePre,FilePost,FileExt,SeqFormat):
	TempDict = defaultdict(dict)#TempDict[Locus][ContigName] = ContigSeq
	NumFiles = 0
	#and make a new locus list for the loci we actually have
	LListNew = [ ]
	for Locus in LList:
		FileName = Folder + FilePre + Locus + FilePost + FileExt
		try:
			InFile = open(FileName, 'rU')
			for record in SeqIO.parse(InFile, SeqFormat):
				TempDict[Locus][record.id] = str(record.seq)
			#print("%d sequences for locus %s were read from the file %s.\n" % (len(TempDict[Locus].keys()), Locus, FileName))
			InFile.close()
			NumFiles += 1
			LListNew.append(Locus)
		except IOError:
			"do nothing"
			#print("There is no file for locus %s.\n" % (Locus))
	print("%d sequence files were read, with names such as %s.\n" % (NumFiles, FileName))
	#sys.stderr.write("%d sequence files were read, with names such as %s.\n" % (NumFiles, FileName))
	return TempDict, LListNew
	#TempDict is ContigDict.
	#LListNew is LocusList (after LocusSeqGetter is run).

#TreeNumReader reads a guide tree output by RAxML with [I##] node and tip numbers.
#It is called something starting with RAxML_originalLabelledTree
#The tip numbers are added to a dictionary (DictTemp), while the node numbers
#are kept as the node numbers on the tree.
def TreeNumReader(TFName,TreeFormat):
	#various regular expressions that are needed to read this tree format
	MyRe1 = r"\):(\d+\.\d+)\[I(\d+)\]"
	MySub1 = r")\2:\1"
	MyRe2 = r":(\d+\.\d+)\[I\d+\]"
	MySub2 = r":\1"
	MyRe3 = r"([\w\.-]+):\d+\.\d+\[I(\d+)\]"
	DictTemp = { } #DictTemp[TaxNum] = Taxon
	#reformatting the newick string
	InFile = open(TFName, 'rU')
	for Line in InFile:
		TreeString1 = Line.strip('\n').strip('\r')
	InFile.close()
	TreeString2 = re.sub(MyRe1, MySub1, TreeString1)
	TreeString3 = re.sub(MyRe2, MySub2, TreeString2)
	#and adding the taxon numbers to the dictionary
	TreeSplit = TreeString2.split(',')
	for TreePart in TreeSplit:
		MyResult = re.search(MyRe3, TreePart)
		DictTemp[MyResult.group(2)] = MyResult.group(1)
	#opening the reformatted tree
	tree1 = dendropy.Tree.get_from_string(TreeString3, schema=TreeFormat, preserve_underscores=True, as_rooted=False)
	#The next line prints the tree and can be uncommented, if you want to see it.
	#print(tree1.as_ascii_plot(show_internal_node_labels=True))
	#print len([str(node2.taxon) for node2 in tree1.leaf_nodes()])
	return(DictTemp, tree1)
	#and OTUNumDict and Tree1, respectively, in PosFinder.

#NodeCompare determines the sets of nodes where each contig could be
def NodeCompare(BackTree, NumDict, FileName, LCDict, COVal):
	DictTemp = defaultdict(list) #DictTemp[SeqName] = list of sister nodes		
	SeqName = ""
	InFile = open(FileName, 'rU')
	#reading the RAxML_classificationLikelihoodWeights file and getting the information of which nodes each
	#contig could be sister to
	#for the TestP confidence interval where this contig could be
	SeqListTemp = [ ]
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split(' ')
		SeqName = Line[0]
		SeqListTemp.append(SeqName)
		NodeNum = Line[1][1:]
		PVal = float(Line[2])
		#only add information about this placement to the dictionary if the placement has a high enough probability
		if PVal >= COVal:
			DictTemp[SeqName].append(NodeNum)
	SeqListTemp = list(set(SeqListTemp))
	#if len(SeqListTemp) != len(DictTemp):
		#print("Only %d of the %d contigs had positions supported by probabilities greater than the cutoff value of %.3f.  The remaining %d contigs will not be used.\n" % (len(DictTemp), len(SeqListTemp), COVal, len(SeqListTemp)-len(DictTemp)))
		#sys.stderr.write("Only %d of the %d contigs had positions supported by probabilities greater than the cutoff value of %.3f.  The remaining %d contigs will not be used.\n" % (len(DictTemp), len(SeqListTemp), COVal, len(SeqListTemp)-len(DictTemp)))
	GDictTemp = defaultdict(list) #GDictTemp[GroupNum] = list of Contigs
	NDictTemp = { } #NDictTemp[NodeNum] = GroupNum
	NListDict = defaultdict(list) #NListDict[GroupNum] = list of nodes
	GroupNum = 1
	for SeqName in DictTemp:
		SeqLen = len(LCDict[SeqName])
		#print ("%s: length: %d, groups: %s\n" % (SeqName, SeqLen, ", ".join(DictTemp[SeqName])))
		GroupListTemp = [ ]
		for NodeNum in DictTemp[SeqName]:
			try:
				GroupListTemp.append(NDictTemp[NodeNum])
			except KeyError:
				"do nothing"
		GroupListTemp = list(set(GroupListTemp))
		#if only one group is represented in the nodes:
		if len(GroupListTemp) == 1:
			#print("Only group %d present." % (GroupListTemp[0]))
			#add this sequence to an existing group
			NewGroup = GroupListTemp[0]
			GDictTemp[NewGroup].append(SeqName)
			for NodeNum in DictTemp[SeqName]:
				#and give that group to all remaining nodes
				NDictTemp[NodeNum] = NewGroup
				#and add its nodes to the NListDict:
				NListDict[NewGroup].append(NodeNum)
		#but if none of the nodes have a group yet
		elif len(GroupListTemp) == 0:
			#assign them the next group
			#print("No groups found for this contig.  It will become Group %d." % (GroupNum))
			GDictTemp[GroupNum] = [SeqName]
			for NodeNum in DictTemp[SeqName]:
				NDictTemp[NodeNum] = GroupNum
				NListDict[GroupNum].append(NodeNum)
			GroupNum += 1
		else:
			#print("This contig had several groups: %s.  They will be combined under group %d." % (",".join([str(GN) for GN in GroupListTemp]), GroupListTemp[0])) 
			#arbitrarily choose the first group number
			NewGroup = GroupListTemp[0]
			NListDict[NewGroup] += DictTemp[SeqName]
			GDictTemp[NewGroup].append(SeqName)
			#combine the existing groups in GDictTemp
			for GNumTemp in GroupListTemp:
				if GNumTemp != NewGroup:
					GDictTemp[NewGroup] += GDictTemp[GNumTemp]
					#get rid of the extra entries in GDictTemp
					del(GDictTemp[GNumTemp])
					#do the same thing for the dictionary that says which nodes belong to which groups
					NListDict[NewGroup] += NListDict[GNumTemp]
					del(NListDict[GNumTemp])
			NListDict[NewGroup] = list(set(NListDict[NewGroup]))
			#add this group to all nodes
			for NodeNum in NListDict[NewGroup]:
				NDictTemp[NodeNum] = NewGroup
	for GroupNum in NListDict:
		NListDict[GroupNum] = list(set(NListDict[GroupNum]))
	return (GDictTemp, NListDict)
	#DictTemp could become ContigPosDict, which says which nodes a sequences could be sister to.
	#GDictTemp could become GroupDict, which says which contigs are members of which groups.
	#NListDict could become NodeDict, which says which nodes are classified in which groups.

#GroupMerge recircumscribes the groups found in GDict if they overlap in CDict
def GroupMerging(GDictIn, CDictIn):
	#GDict[GroupNum] = list of group members
	#CDictIn[GroupNum] = list of items in group--do not want groups to have overlapping items
	CtoGDict = defaultdict(list) #CtoGDict[Item] = list of GroupNums
	GCombineDict = defaultdict(list) #GCombineDict[GroupNum] = list of other GroupNums that overlap with it
	NGDict = { } #NGDict[NGroupNum] = list of old GroupNums
	OGDict = { } #OGDict[GroupNum] = NGroupNum it is part of
	GDictOut = defaultdict(list) #GDictOut[NGroupNum] = list of group members
	CDictOut = defaultdict(list) #CDictOut[NGroupNum] = list of items
	#reversing CDictIn
	for GroupNum in CDictIn:
		for Item in CDictIn[GroupNum]:
			CtoGDict[Item].append(GroupNum)
	#now looking through the reversed dictionary to see if any of the items have more than one group
	for Item in CtoGDict:
		#add the groups to GCombineDict for each GroupNum
		for GroupNum in CtoGDict[Item]:
			GCombineDict[GroupNum] += CtoGDict[Item]
	#now going through that dictionary and condensing it, so each group is only listed once
	for GroupNum in GCombineDict:
		GCombineDict[GroupNum] = list(set(GCombineDict[GroupNum]))
	for GroupNum in GCombineDict:
		#if this group doesn't overlap with any other group
		if len(GCombineDict[GroupNum]) == 1:
			NGDict[GroupNum] = [GroupNum]
			OGDict[GroupNum] = GroupNum
		else:
			ListTemp = [ ]
			#look for existing group names
			for GN2 in GCombineDict[GroupNum]:
				try:
					ListTemp.append(OGDict[GN2])
				except KeyError:
					"do nothing"
			#condensing the list so that each member is present only once
			ListTemp = list(set(ListTemp))
			#if there is no existing name, make a new group
			if len(ListTemp) == 0:
				NGroupNumTemp = GroupNum
				NGDict[NGroupNumTemp] = GCombineDict[GroupNum]
			#if there is one existing name, use it
			elif len(ListTemp) == 1:
				NGroupNumTemp = ListTemp[0] 
				NGDict[NGroupNumTemp] += GCombineDict[GroupNum]
				NGDict[NGroupNumTemp] = list(set(NGDict[NGroupNumTemp]))
			#if there are multiple names, they need to be combined :(
			elif len(ListTemp) > 1:
				#use the smaller group number
				NGroupNumTemp = sorted(ListTemp)[0]
				#make a list of the members of this combined group
				ListTemp2 = [ ]
				#first, add the group members from the original GCombineDict
				ListTemp2 += GCombineDict[GroupNum]
				#then add all of the members from the old groups that overlapped with this list
				for NGN2 in ListTemp:
					ListTemp2 += NGDict[NGN2]
				#condense the list
				ListTemp2 = list(set(ListTemp2))
				#This list is the new set of members of the combined group
				NGDict[NGroupNumTemp] = ListTemp2
				#now delete the other overlapping groups
				for NGN2 in ListTemp:
					if NGN2 != NGroupNumTemp:
						del NGDict[NGN2]
				#and update OGDict with the new placement of each member
				for NGN2 in ListTemp2:
					OGDict[NGN2] = NGroupNumTemp
			for GN2 in GCombineDict[GroupNum]:
				OGDict[GN2] = NGroupNumTemp
		#then I need to make the new group list
	for NGroupNum in NGDict:
		#if len(NGDict[NGroupNum]) > 1:
			#print("The old groups %s were merged to form the new group %d.\n" % (", ".join(str(GroupNum) for GroupNum in NGDict[NGroupNum]), NGroupNum))
		for GroupNum in NGDict[NGroupNum]:
			GDictOut[NGroupNum] += GDictIn[GroupNum]
	for GroupNum in CDictIn:
		NGroupNum = OGDict[GroupNum]
		CDictOut[NGroupNum] += CDictIn[GroupNum]
	return(GDictOut, CDictOut)	
		
#PosFinder goes through all of the contigs from all of the loci, divides them into groups based on their
#position in the tree, and determines which paralog(s) each group belongs to.
def PosFinder(LList, Folder, FilePre, PDict, CDict):
	CParDict = defaultdict(dict)#CPDict[Locus][Paralog] = List of Contigs
	for Locus in LList:
		CParDict[Locus] = defaultdict(list) #CPDict[Locus][GroupNum]=list of Contigs
		LPDict = defaultdict(list)#LPDict[ParName]=list of groups having that name
		GNameDict = { } #GNameDict[GroupNum] = GroupName
		GoodContigs = 0
		BadContigs = 0
		#read the tree
		#read the tree
		TreeFileName = Folder+"RAxML_originalLabelledTree."+FilePre+Locus
		(OTUNumDict, Tree1) = TreeNumReader(TreeFileName, 'newick')
		#read the likelihood file and determine which nodes a sequence could potentially be located at
		LikeFileName = Folder+"RAxML_classificationLikelihoodWeights."+FilePre+Locus
		(GroupDict, NumGroupDict) = NodeCompare(Tree1, OTUNumDict, LikeFileName, CDict[Locus], CutOff)
		#print("The contigs for locus %s were in %d groups.\n" % (Locus, len(NumGroupDict)))
		SisterGroupDict = defaultdict(list)
		for GroupNum in NumGroupDict:
			#print("Group %d contained %d contigs from the following nodes: %s.\n" % (GroupNum, len(GroupDict[GroupNum]), ", ".join(NumGroupDict[GroupNum])))
			for NodeNum in NumGroupDict[GroupNum]:
				try:
					SisterGroupDict[GroupNum].append(OTUNumDict[NodeNum])
					#print ("%s: %s\n" % (NodeNum, OTUNumDict[NodeNum]))
				except KeyError:
					node = Tree1.find_node_with_label(NodeNum)
					ListTemp = [str(node2.taxon) for node2 in node.leaf_nodes()]
					SisterGroupDict[GroupNum]+= ListTemp
					#print ("%s: %s" % (NodeNum, ", ".join(ListTemp)))
			SisterGroupDict[GroupNum] = list(set(SisterGroupDict[GroupNum]))
		#merging overlapping groups with overlapping sister groups, if desired
		if GrpFate == "merge":
			(GroupDict, SisterGroupDict) = GroupMerging(GroupDict, SisterGroupDict)
			#print("After merging, the contigs for locus %s were in %d groups.\n" % (Locus, len(GroupDict)))
		#changing the names if two groups belong to the same paralog
		for GroupNum in SisterGroupDict:
			ParList = [ ]
			for Sister in SisterGroupDict[GroupNum]:
				try:
					ParList.append(PDict[Locus][Sister])
					#print("%s: %s\n" % (Sister, PDict[Locus][Sister]))
				except KeyError:
					"no paralogs for this sequence"
					#print("%s: no paralogs listed\n" % (Sister))
			ParList = list(set(ParList))
			#Determining if a group of sequences can be classified according to paralog and naming it accordingly.
			if len(ParList) == 0:
				#print("None of the sister sequences were classified as to paralog.\n")
				ParNameTemp = "Ambig_"+Locus+"_none"
			elif len(ParList) == 1:
				#print("This group belongs to the paralog %s.\n" % (ParList[0]))
				ParNameTemp = ParList[0]
			elif len(ParList) == 2:
				#print("This group could belong to either of two paralogs: %s or %s.\n" % (ParList[0], ParList[1]))
				ParNameTemp = "Ambig_"+Locus+"_"+"_".join(ParList)
			elif len(ParList) > 2:
				#print("It is not clear which of the following paralogs this group belongs to: %s.\n" % (", ".join(ParList)))
				ParNameTemp = "Ambig_"+Locus+"_"+str(len(ParList))
			LPDict[ParNameTemp].append(GroupNum)
			GNameDict[GroupNum] = ParNameTemp
			#if len(GroupDict[GroupNum]) == 1:
				#print("The contig is %s, with a length of %d.\n" % (GroupDict[GroupNum][0], len(CDict[Locus][GroupDict[GroupNum][0]])))
		for ParNameTemp in sorted(LPDict.keys()):
			#print("%d groups were found for paralog %s.\n" % (len(LPDict[ParNameTemp]), ParNameTemp))
			if len(LPDict[ParNameTemp]) > 1:
				#Then ParNameTemp is not enough of a name for those groups.  They need to get a number as well.
				Num = 1
				#looking through all the groups
				for GroupNum in GNameDict:
					#if the name of the group is ParNameTemp,
					if GNameDict[GroupNum] == ParNameTemp:
						#giving it a number as well
						NameTemp = GNameDict[GroupNum]+"_P"+str(Num)
						GNameDict[GroupNum] = NameTemp
						#and increasing the number by one for the next renaming
						Num += 1
		#adding the contigs for that locus to the CParDict:
		for GroupNum in GroupDict:
			LocusGroup = GNameDict[GroupNum]
			CParDict[Locus][LocusGroup] = GroupDict[GroupNum]
	return CParDict
	#This is ContigParDict

#ParalogSeqWriter writes sequences in the desired format to files from two dictionaries.  The first
#dictionary is a dictionary of the sequences classified by locus (and contig name), while the second
#dictionary does not have the sequences, but says which contigs belong to each locus.
#modified from LocusParalogSeqWriter from tclade_finder.py,
#which was modified from LocusGroupSeqWriter from tbaits_intron_removal.py
def ParalogSeqWriter(CDict, CPDict, Folder, Prefix, Suffix, SeqFormat):
	#CPDict[Locus][Paralog]=list of Contigs
	#CDict[Locus][ContigName] = ContigSeq
	TempDict = defaultdict(dict)#TempDict[Locus][Paralog] = file name
	for Locus in CPDict:
		NumFiles = 0
		for Paralog in CPDict[Locus]:
			OutFileName = Folder+Prefix+Paralog+Suffix
			TempDict[Locus][Paralog] = Prefix+Paralog
			OutFile = open(OutFileName, 'w')
			for ContigName in CPDict[Locus][Paralog]:
				Record1 = SeqRecord(seq=Seq(CDict[Locus][ContigName], IUPAC), id = ContigName, description = "")
				SeqIO.write(Record1, OutFile, SeqFormat)
			OutFile.close()
			NumFiles += 1
		print("%d sequence files were written for the locus %s, with names such as %s.\n" % (NumFiles, Locus, OutFileName))
		#sys.stderr.write("%d sequence files were written for the locus %s, with names such as %s.\n" % (NumFiles, Locus, OutFileName))
	return TempDict
	#This is FileDict
	
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
			Line += "mafft --addfragments "+Folder+NamePart+".fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+NamePart+"_pars_al.fa\n"
			OutList.append(Line)
	OutFileName = Folder+Prefix+"Analysis_Script.sh"
	OutFileWriting(OutFileName, OutList)
	print("The shell script for analyzing the sequences further was written to %s.\n" % (OutFileName))
	sys.stderr.write("The shell script for analyzing the sequences further was written to %s.\n" % (OutFileName))
	
#read the names of the sequence files from blast file list read in tbaits_intron_removal.py
LocusList = CaptureColumn(LocusFileList, 0)

#make a dictionary of which backbone sequences belong to which paralogs
if PFileName != "none":
	ParalogDict = DDictFromFile(PFileName)
else:
	print("No file was given for classifying the backbone sequences into paralogs.  It is assumed that their classification is unknown.\n")
	sys.stderr.write("No file was given for classifying the backbone sequences into paralogs.  It is assumed that their classification is unknown.\n")

#read the sequence files
(ContigDict, LocusList) = LocusSeqGetter(LocusList, SeqFolder, SeqFilePre, SeqFilePost, '.fa', 'fasta')

#determining which paralogs the contigs belong to
ContigParDict = PosFinder(LocusList, TreeFolder, TreeFilePre, ParalogDict, ContigDict)

#write the output
#sequence files for each group/locus:
FileDict = ParalogSeqWriter(ContigDict, ContigParDict, OutFolder, OutFilePre, '.fa', 'fasta')

#script for analyzing the sequences
AMScriptWriter(FileDict, OutFolder, OutFilePre, AlFolder, AlFilePre, AlFilePost)

#an information file:
#Locus\tContig\tGroupNum\tParalog
OutGroupList = ['Locus\tContig\tLength\tParalog\tLength\n']
for Locus in ContigParDict:
	for Paralog in ContigParDict[Locus]:
		for Contig in ContigParDict[Locus][Paralog]:
			LineG = Locus+"\t"+Contig+"\t"+Paralog+"\t"+str(len(ContigDict[Locus][Contig]))+"\n"
			OutGroupList.append(LineG)
OutFileGName = OutFolder+OutFilePre+"Contig_Groups.txt"
OutFileWriting(OutFileGName, OutGroupList)
