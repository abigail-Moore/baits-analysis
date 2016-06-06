#! /usr/bin/env python

#tcontig_selection.py version 1.1 17 Feb. 2016 Abby Moore
#This script follows tcontigs_to_fixed_paralogs.py and goes through the alignments of contigs that can't be 
#assembled into good enough consensus sequences.  It selects the best segment of the alignment using four criteria
#involving the contigs with significant overlap with that region:
#number of contigs, number of individuals, mean length of overlap, median total contig length.
#The contigs that have significant overlap with the chosen segment are then used in their entirety, not just in the
#overlapping region.
#It reads the list of paralogs and groups from the OutFilePre_Contigs_to_Redo.txt file, produced by tundivcontigs_combiner.py
#script by combining the undivisible contigs from the various groups.
#Version 1.1 is modified so that it only selects contigs over 100bp.

'''
(no header, tab-delimitted)
(one line per group, even when everything else is the same)
amk: Locus [0]
amk1: Paralog [1]
Montiaceae: Group[2]
list of contig names (comma-separated): ContigNameList [3]
'''
#It also needs a Groups_List file:
'''
(no header, tab-delimitted)
Anacampserotaceae: Group [0]
Anacs/none: GroupName (or "none") [1]
'''

'''
tcontig_selection.py PFileName GFileName InFolder InFilePre OutFolder OutFilePre SeqFilePath SeqFolderPre SeqFolderPost SeqFilePre IndDictFileName
/gpfs/scratch/ajm3/eedwards/scripts/tcontig_selection.py /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_contigsplit/Ln1tcs_Contigs_to_Redo.txt /gpfs/scratch/ajm3/eedwards/general/Group_List_round1_all.txt /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_contigsplit/ Ln1tcs_ same Ln1tcs_ /gpfs/scratch/ajm3/eedwards/ s2_ spades_Ln1t_sequences/ Ln1ts5_ /gpfs/scratch/ajm3/eedwards/general/Ln1_inds2.txt
'''

import sys
import numpy #to calculate means and medians
from Bio import AlignIO #to read sequence alignments
from Bio.Align import MultipleSeqAlignment #to make sequence alignments
from collections import defaultdict #to make dictionaries with multiple levels
import random #to choose items at random from a list

Usage = '''
tcontig_selection.py version 1.1
This script looks through alignments with multiple contigs per individual and the
name of the individual separated from the rest of the name of the contig by a "-"
symbol.  It chooses a subset of contigs to analyze further for each individual,
such that the contigs all overlap in the alignment.  When multiple such sets of 
contigs exist, it takes the set of the longest contigs.
[file with information about the paralogs: Contigs_to_Redo.txt file from
tundivcontigs_combiner.py]
[list of groups--tab delimitted, group (tab) abbreviated group name]
[folder in which alignments are found]
[prefix for alignment files--everything before the paralog name]
[folder to which to output the files with all sequences for each locus, or 
"same", if the same as the input folder]
[prefix for the combined output files--everything before the paralog name]
[directory in which the folders for the individual output files are found]
[anything between the name of the main folder and the name of the subfolder or 
"none", if nothing]
[anything after the name of the subfolder, or "none", if nothing--this can just
be "/", if there is no suffix, but there is still a folder]
[prefix for the sequence files, or "none" if none]
[dictionary saying which group each individual belongs to]
'''

if len(sys.argv) != 12:
	sys.exit("ERROR!  tcontig_selection.py expects 11 addition arguments and you supplied %d!\n%s" % (len(sys.argv)-1, Usage))
PFileName = sys.argv[1]
GFileName = sys.argv[2]
InFolder = sys.argv[3]
if InFolder[-1] != "/":
	InFolder += "/"
InFilePre = sys.argv[4]
OutFolder = sys.argv[5]
if OutFolder == "same":
	OutFolder = InFolder
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[6]
SeqFilePath = sys.argv[7]
if SeqFilePath[-1] != "/":
	SeqFilePath += "/"
SeqFolderPre = sys.argv[8]
if SeqFolderPre == "none":
	SeqFolderPre = ""
SeqFolderPost = sys.argv[9]
if SeqFolderPost[-1] != "/":
	SeqFolderPost += "/"
if SeqFolderPost == "none/":
	SeqFolderPost = ""
SeqFilePre = sys.argv[10]
if SeqFilePre == "none":
	SeqFilePre = ""
IndDictFileName = sys.argv[11]

Verbose = False

#DictFromFile makes a dictionary from a tab-delimited file, where the keys are in columns KC, and
#the values are in column VC
#a modified version of the function in from tbaits_intron_removal.py
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
	#This is GroupDict and IndDict

#DListDictFromFile makes a four-level dictionary from a tab-delimited file, where the last level is a list
#and the columns to go in each level can be specified
#modified from DDictFromFile in tseq_placer_dup.py
def DListDictFromFile(FileName,Key1Col,Key2Col,Key3Col,ListCol,Delim):
	TempDict = defaultdict(dict)
	LineNum = 0
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		ListTemp = Line[ListCol].split(Delim)
		try:
			TempDict[Line[Key1Col]][Line[Key2Col]][Line[Key3Col]] = ListTemp
		except KeyError:
			TempDict[Line[Key1Col]][Line[Key2Col]] = defaultdict(dict)
			TempDict[Line[Key1Col]][Line[Key2Col]][Line[Key3Col]] = ListTemp
		LineNum += 1
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\n" % (LineNum, FileName))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\n" % (LineNum, FileName))
	return TempDict
	#This is ParalogDict

#AlignmentSubsetting takes an alignment and makes a new alignment that only has the sequences in a given list
def AlignmentSubsetting(AlignOrig, SeqList):
	#making a new alignment with just our contigs
	PrunedAlign = ""
	for record in AlignOrig:
		#if it is one of our contigs
		if record.id in SeqList:
			#if so, try to add that contig to the new alignment
			try:
				PrunedAlign.append(record)
			#but, if this is the first contig we want, we need to initialize the new alignment with that contig instead
			except AttributeError:
				PrunedAlign = MultipleSeqAlignment([record])
	return PrunedAlign
	#This is PAlignGaps for the whole script and other thing in functions.

#GapFinder finds the gaps in an alignment.  It returns a dictionary in which the gaps are numbered consecutively
#(starting with 0) of the following form: DictTemp[GapNum]['start'/'end']: position in the alignment where the gap starts
#or ends, respectively
def GapFinder(AlignTemp):
	#look through each position of the alignment
	DictTemp = defaultdict(dict)#DictTemp[GapNum]['start'/'end'] = SeqPos
	GapNum = 0
	GapStart = 0
	InGap = True
	GapLength = 0
	for SeqPos in range(0, AlignTemp.get_alignment_length()):
		#determine whether or not there are nucleotides at that position
		GapOnly = True
		for record in AlignTemp:
			if record[SeqPos] != '-':
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
		#in the middle of a segment [(GapOnly == False) and (InGap == False)]
		#or in the middle of a gap [(GapOnly == True) and (InGap == True)]
	#check to see if there is a gap at the end of the alignment
	if (InGap == True):
		DictTemp[GapNum]['start'] = GapStart
		DictTemp[GapNum]['end'] = SeqPos
		GapLength += (SeqPos-GapStart+1)
	if Verbose == True: print("%d gaps were found, with a total length of %d.\n" % (GapNum+1, GapLength))
	return DictTemp
	#This is GapDict in AlignmentGapRemoving

#AlignmentGapRemoving takes an alignment and removes all of the characters that are gaps
#in all sequences.  It calls GapFinder.
def AlignmentGapRemoving(AlignGaps):
	GapDict = GapFinder(AlignGaps)
	AlEnd = AlignGaps.get_alignment_length()-1
	if Verbose == True: print("Before removing any gaps, the alignment is %d bases long.\n" % (AlEnd+1))
	for GapNum in sorted(GapDict.keys(), reverse=True):
		#print("Before removing gap %d: gap length: %d, alignment length %d\n" % (GapNum, GapDict[GapNum]['end']-GapDict[GapNum]['start']+1, AlignGaps.get_alignment_length()))
		#if this gap encompasses the end of the alignment, cut that bit off of the end of the alignment
		if GapDict[GapNum]['end'] == AlEnd:
			AlignGaps = AlignGaps[:, :GapDict[GapNum]['start']]
		#if this gap encompasses the start of the alignment, cut off the start
		elif GapDict[GapNum]['start'] == 0:
			AlignGaps = AlignGaps[:, (GapDict[GapNum]['end']+1):]
		#for all gaps in the middle of the alignment
		else:
			AlignGaps = AlignGaps[:, :GapDict[GapNum]['start']] + AlignGaps[:, (GapDict[GapNum]['end']+1):]
		#print("After removing gap %d: alignment length %d\n" % (GapNum, AlignGaps.get_alignment_length()))
	if Verbose == True: print("After removing %d gaps, the alignment is now %d bases long.\n" % (len(GapDict), AlignGaps.get_alignment_length()))
	return AlignGaps
	#This is now PAlign

#AlInfoFinder makes three things:
#a dictionary containing information about sequences in an alignment (DictTemp)
#a dictionary with the list of contigs that belong to each individual (INDict)
#a dictionary with the list of contigs that overlap every other contig (OVDict)
def AlInfoFinder(AlignTemp, SepChar, Overlap):
	DictTemp = defaultdict(dict)#DictTemp[Contig]['IndName'/'CStart'/'CEnd'/'CLength']
	INDict = defaultdict(list)#INDict[IndName] = list of Contigs
	OVDict = defaultdict(list)#OVDict[Contig] = list of Contigs with which it overlaps
	AlLength = AlignTemp.get_alignment_length()
	AlRange = range(0, AlLength)
	for record in AlignTemp:
		Contig = record.id
		IndName = record.id.split(SepChar)[0]
		DictTemp[Contig]['IndName'] = IndName
		INDict[IndName].append(Contig)
		CStart = -1
		CEnd = -1
		CLength = 0
		#go through the sequence
		for SeqPos in AlRange:
			#the base is not a gap,
			if record[SeqPos] != "-":
				#increase the sequence length by 1
				CLength += 1
				#assume that the sequence also ends here
				CEnd = SeqPos
				#and if this is the first non-gap base, it is the start of the sequence, as well
				if CStart == -1:
					CStart = SeqPos
		#record the information to DictTemp
		DictTemp[Contig]['CStart'] = CStart
		DictTemp[Contig]['CEnd'] = CEnd
		DictTemp[Contig]['CLength'] = CLength
		Length2 = CEnd-CStart+1
		if Verbose == True:
			if Length2 > 1.25*CLength:
				print("WARNING! Contig %s has many gaps.  It covers %d bases in the alignment, but only has information at %d of those bases.\n" % (Contig, Length2, CLength))
			print("%s: start: %d, end: %d, length1: %d, length2: %d\n" % (Contig, CStart, CEnd, CLength, Length2))
	#now looking through the contigs again, to determine overlap
	for Contig in DictTemp:
		CS1 = DictTemp[Contig]['CStart']
		CE1 = DictTemp[Contig]['CEnd']
		CL1 = DictTemp[Contig]['CLength']
		for Contig2 in DictTemp:
			CS2 = DictTemp[Contig2]['CStart']
			CE2 = DictTemp[Contig2]['CEnd']
			CL2 = DictTemp[Contig2]['CLength']
			OVS = max(CS1, CS2)
			OVE = min(CE1, CE2)
			OVL = OVE-OVS+1
			if (OVL > Overlap*CL1) or (OVL > Overlap*CL2):
				if Verbose == True: print("Sequences %s (length %d) and %s (length %d) overlap by %d bases.\n" % (Contig, CL1, Contig2, CL2, OVL))
				OVDict[Contig].append(Contig2)
	return DictTemp, INDict, OVDict
	#DictTemp is CInfoDict
	#INDict is IndSeqDict
	#OVDict is OverDict

#BreakFinder analyzes the results of AlInfoFinder (specifically DictTemp) to find
#places in the alignment where there are unusually many starts and ends of sequences
def BreakFinder(IDict, IntPercent):
	SEDict = { } #the dictionary of contig starts and ends
	for Contig in IDict:
		CS = IDict[Contig]['CStart']
		CE = IDict[Contig]['CEnd']
		try: SEDict[CS] += 1
		except KeyError: SEDict[CS] = 1 
		try: SEDict[CE] += 1
		except KeyError: SEDict[CE] = 1
	if Verbose == True: print("%s" % (", ".join([str(Num)+": "+str(SEDict[Num]) for Num in sorted(SEDict.keys())])))
	SEBPDict = defaultdict(dict)
	AlEnd = sorted(SEDict.keys())[-1]
	IntervalSize = IntPercent*AlEnd
	if IntervalSize > 98:
		IntervalSize = 98
	if Verbose == True: print("Starts and ends of contigs within %.2f of one another, will be considered to be starting or ending at the same breakpoint.\n" % (IntervalSize))
	#There is also the question if I want more information in the SEDict, such as a list of contigs (the dictionary could even be the list of contigs, and
	#the number of individuals with each break point could be the length of the list--this would eliminate the KeyError stuff, although perhaps I could use defaultdict(int) and
	#get rid of it that way instead).
	BPDict = defaultdict(dict)
	BPStart = -99
	BPEnd = -99
	BPNum = 0
	for SE in sorted(SEDict.keys()):
		if SE > BPStart + IntervalSize:
			BPStart = SE
			BPNum += 1
			BPDict[BPNum]['BPStart'] = SE
			BPDict[BPNum]['BPMems'] = SEDict[SE]
			BPDict[BPNum]['BPEnd'] = SE
			BPDict[BPNum]['BPList'] = [SE]
		else:
			BPDict[BPNum]['BPEnd'] = SE
			BPDict[BPNum]['BPMems'] += SEDict[SE]
			BPDict[BPNum]['BPList'].append(SE)
	#Deciding which breakpoints involve enough of the sequences to be useful
	BPMemList = [ ]
	for BPNum in BPDict:
		BPMemList.append(BPDict[BPNum]['BPMems'])
	#The question here is whether we want to use the mean of the list of BPMems or the mean of the list(set(BPMems)).
	#The latter eliminates all repeated BPMems numbers, which basically shifts the mean to the right, because the
	#repeated BPMems numbers are mainly 1s, 2s, and 3s.
	BPThreshold = numpy.median(list(set(BPMemList)))
	if Verbose == True: print("The alignment will be broken up at all break points where at least %.1f contigs start or end.\n" % (BPThreshold))
	#Finding the positions of all of the breakpoints that are common enough.
	MajBPList = [ ]
	for BPNum in BPDict:
		if BPDict[BPNum]['BPMems'] >= BPThreshold:
			BPPair = [BPDict[BPNum]['BPStart'], BPDict[BPNum]['BPEnd']]
			MajBPList.append(BPPair)
	if Verbose == True: print("There are %d such breakpoints.\n" % (len(MajBPList)))
	MajBPList = sorted(MajBPList)
	#Finding the segments of the sequence
	DictTemp = defaultdict(dict)#DictTemp[SegNum]['SegStart'/'SegEnd'] = SeqPos (this will be built upon later)
	SegNum = 1
	#if the contigs start in too many different places for the start of the alignment to count as the first breakpoint,
	#including those anyway
	if MajBPList[0][0] > 10:
		DictTemp[SegNum]['SegStart'] = 0
		SegNum += 1
	#Now going through the list of break points and adding those to the dictionary
	#for each break point range,
	for BPPair in MajBPList:
		#the end of the range is the start of the current alignment segment
		DictTemp[SegNum]['SegStart'] = BPPair[1]
		#and the start of the range is the end of the previous alignment segment
		if SegNum > 1:
			DictTemp[SegNum-1]['SegEnd'] = BPPair[0]
		SegNum += 1
	#now determining if the last fully-written ended at the end of the alignment:
	#if not, making a breakpoint that does
	if DictTemp[SegNum-1]['SegStart'] < (AlEnd-10):
		DictTemp[SegNum-1]['SegEnd'] = AlEnd
	#if so, deleting the partially-written breakpoint that would have gone over the edge of the alignment
	else:
		del DictTemp[SegNum-1]
	if Verbose == True: print("These breakpoints divided the alignment into %d segments.\n" % (len(DictTemp)))
	return DictTemp
	#This is AlSegDict

#SegMembershipFinder takes a dictionary of alignment segments and a dictionary showing which contigs are where in the
#alignment, and determines which contigs and which individuals are in each alignment segment.
#This works one segment per alignment (for example one particular amino acid of interest) or multiple segments (for example,
#trying to find the highest-coverage segment in the alignment).
def SegMembershipFinder(SegDictIn, CDict, OVPercent, PGDict):
	#CDict[Contig]['IndName'/'CStart'/'CEnd'/'CLength']
	#before: SegDictIn[SegNum]['SegStart'/'SegEnd'] = SeqPos
	#after: SegDictIn[SegNum]['SegStart'/'SegEnd'/'ContigList'/'NumContigs'/'NumInds'/'MeanOVLength'/'MeanCLength'/'TotOVLen'/
	#'IndDict'([IndName] = NumContigs)/'Overlap'(list with amount of overlap of each contig)/
	#'ConLength'(list with the total sequence length of each contig)]
	ISDict = defaultdict(int) #ISDict[IndName] = TotalNumSegments
	for SegNum in SegDictIn:
		SegDictIn[SegNum]['IndDict'] = defaultdict(int) #SegDictIn[SegNum]['IndDict'][IndName] = NumContigs
		SegDictIn[SegNum]['ContigList'] = [ ]#list with the contigs that include this region
		SegDictIn[SegNum]['Overlap'] = [ ]#list with the amounts each contig overlaps this region
		SegDictIn[SegNum]['ConLength'] = [ ]#list with the lengths of each of the contigs (not just the overlapping bits)
		SegDictIn[SegNum]['GroupList'] = [ ]#list with the groups that have that segment
		SegStart = SegDictIn[SegNum]['SegStart']
		SegEnd = SegDictIn[SegNum]['SegEnd']
		SegLength = SegEnd - SegStart + 1
		for Contig in CDict:
			OVStart = max(SegStart, CDict[Contig]['CStart'])
			OVEnd = min(SegEnd, CDict[Contig]['CEnd'])
			OVLength = OVEnd - OVStart + 1
			#if the overlap is long enough, accept the contig, and add it to the dictionary
			if OVLength > SegLength*OVPercent:
				IndName = CDict[Contig]['IndName']
				GroupName = PGDict[Contig]
				SegDictIn[SegNum]['IndDict'][IndName] += 1
				SegDictIn[SegNum]['ContigList'].append(Contig)
				SegDictIn[SegNum]['Overlap'].append(OVLength)
				SegDictIn[SegNum]['ConLength'].append(CDict[Contig]['CLength'])
				SegDictIn[SegNum]['GroupList'].append(GroupName)
				ISDict[IndName] += 1 
	#checking to make sure that nothing unexpected happened, and calculating some summary statistics
	for SegNum in SegDictIn:
		NumContigss = len(SegDictIn[SegNum]['ContigList'])
		SegDictIn[SegNum]['NumContigs'] = NumContigss
		NumInds = len(SegDictIn[SegNum]['IndDict'])
		SegDictIn[SegNum]['NumInds'] = NumInds
		ListTemp = SegDictIn[SegNum]['GroupList']
		ListTemp = list(set(ListTemp))
		SegDictIn[SegNum]['GroupList'] = ListTemp
		NumGroups = len(SegDictIn[SegNum]['GroupList'])
		SegDictIn[SegNum]['NumGroups'] = NumGroups
		if NumInds != 0:
			MeanOVLength = numpy.mean(SegDictIn[SegNum]['Overlap'])
			MedCLength = numpy.median(SegDictIn[SegNum]['ConLength'])
		else:
			MeanOVLength = 0
			MedCLength = 0
		SegDictIn[SegNum]['MeanOVLength'] = MeanOVLength
		SegDictIn[SegNum]['MedCLength'] = MedCLength
		TotOVLen = NumContigss*MeanOVLength
		SegDictIn[SegNum]['TotOVLen'] = TotOVLen
		#I am not sure what the most relevant statistic to calculate for the contig length is.
		#Median is quite close to MeanOVLength.  Mean is very biased by a few long contigs.  75th Percentile might be good,
		#but seems a bit arbitrary.
		#Perc75CLength = numpy.percentile(SegDictIn[SegNum]['ConLength'], 75)
		#MeanCLength = numpy.mean(SegDictIn[SegNum]['ConLength'])
		#print("For segment number %d, there were %d overlapping contigs from %d individuals, with a mean overlap length of %.1f, and a median total length of %.1f.\n" % (SegNum, NumContigss, NumInds, MeanOVLength, MedCLength))
		if len(SegDictIn[SegNum]['Overlap']) != NumContigss:
			sys.exit("ERROR!!!!  For segment %d, %d contigs had calculated overlap amounts, but only %d of these are listed in the contig list!!" % (SegNum, len(SegDictIn[SegNum]['Overlap']), NumContigss))
		NumContigsi = 0
		for IndName in SegDictIn[SegNum]['IndDict']:
			NumContigsi += SegDictIn[SegNum]['IndDict'][IndName]
		if NumContigsi != NumContigss:
			sys.exit("ERROR!!!  For segment %d, %d contigs were listed, but only %d of them were classified to individual!!" % (SegNum, NumContigss, NumContigsi))
		if Verbose == True:
			print("%d: SeqStart: %d, SeqEnd: %d, Overlaps: %s, Groups: %d" % (SegNum, SegDictIn[SegNum]['SegStart'], SegDictIn[SegNum]['SegEnd'], ", ".join([str(OVLength) for OVLength in SegDictIn[SegNum]['Overlap']]), SegDictIn[SegNum]['NumGroups']))
			print("%s" % (", ".join([IndName+": "+str(SegDictIn[SegNum]['IndDict'][IndName]) for IndName in sorted(SegDictIn[SegNum]['IndDict'].keys())])))
			print("%s" % (", ".join(sorted(SegDictIn[SegNum]['ContigList']))))
	return SegDictIn, ISDict
	#SegDictIn is SegMemDict.
	#ISDict is IndSegDict.

#SegSelector takes the dictionary with information about the various segments from BreakFinder and SegMembershipFinder and selects the best segment.
#It selects the segment based on four criteria (equally weighted, although this can be changed): the number of contigs that overlap significantly
#with that segment, the number of individuals that have one of those contigs, the mean length of overlap of the overlapping contigs with that segment (so basically
#segment length, assuming the contigs overlap for a lot of the segment), the median (total) length of each overlapping contig, and the total overlap length (number
#of overlapping contigs times mean length of overlap).
def SegSelector(SegDictIn, TotalGroups):
	BestSeg = 0
	#Finding out what the maximum values in each of the five categories is
	MostContigs = 0
	MostInds = 0
	MostOV = 0
	LongestC = 0
	MostTotOV = 0
	for SegNum in sorted(SegDictIn.keys()):
		if SegDictIn[SegNum]['NumContigs'] > MostContigs:
			MostContigs = SegDictIn[SegNum]['NumContigs']
		if SegDictIn[SegNum]['NumInds'] > MostInds:
			MostInds = SegDictIn[SegNum]['NumInds']
		if SegDictIn[SegNum]['MeanOVLength'] > MostOV:
			MostOV = SegDictIn[SegNum]['MeanOVLength']
		if SegDictIn[SegNum]['MedCLength'] > LongestC:
			LongestC = SegDictIn[SegNum]['MedCLength']
		if SegDictIn[SegNum]['TotOVLen'] > MostTotOV:
			MostTotOV = SegDictIn[SegNum]['TotOVLen']
	#Making cutoff values for each category
	#*********These values could potentially be changed.***********
	GSCon = MostContigs - 2
	GSInds = MostInds - 2
	GSOV = MostOV - 30
	GSLongest = LongestC - 30
	GSTotOV = MostTotOV - 150
	#determining whether each segment is above the cutoff value in any of the four categories
	SegScoreDict = defaultdict(list)#SegScoreDict[SegNum] = list of categories in which it was near the best
	AllGroupsList = [ ]
	for SegNum in sorted(SegDictIn.keys()):
		#only looking at segments with a median contig length over 150, because these are the only ones that will be in
		#the final tree
		if SegDictIn[SegNum]['MedCLength'] >= 150:
			#if it has at least the minimum number of overlapping contigs:
			if SegDictIn[SegNum]['NumContigs'] >= GSCon:
				#adding "NumContigs" to the list of categories in which it is near the best
				SegScoreDict[SegNum].append("NumContigs")
			if SegDictIn[SegNum]['NumInds'] >= GSInds:
				SegScoreDict[SegNum].append("NumInds")
			if SegDictIn[SegNum]['MeanOVLength'] >= GSOV:
				SegScoreDict[SegNum].append("MeanOVLength")
			if SegDictIn[SegNum]['MedCLength'] >= GSLongest:
				SegScoreDict[SegNum].append("MedCLength")
			if SegDictIn[SegNum]['TotOVLen'] >= GSTotOV:
				SegScoreDict[SegNum].append("TotOVLen")
			if SegDictIn[SegNum]['NumGroups'] == TotalGroups:
				SegScoreDict[SegNum].append("AllGroups")
				AllGroupsList.append(SegNum)
	#But I guess if there are no segments with median length over 150, then we should run it again to get something.
	if SegScoreDict == { }:
		for SegNum in sorted(SegDictIn.keys()):
			if SegDictIn[SegNum]['NumContigs'] >= GSCon:
				SegScoreDict[SegNum].append("NumContigs")
			if SegDictIn[SegNum]['NumInds'] >= GSInds:
				SegScoreDict[SegNum].append("NumInds")
			if SegDictIn[SegNum]['MeanOVLength'] >= GSOV:
				SegScoreDict[SegNum].append("MeanOVLength")
			if SegDictIn[SegNum]['MedCLength'] >= GSLongest:
				SegScoreDict[SegNum].append("MedCLength")
			if SegDictIn[SegNum]['TotOVLen'] >= GSTotOV:
				SegScoreDict[SegNum].append("TotOVLen")
			if SegDictIn[SegNum]['NumGroups'] == TotalGroups:
				SegScoreDict[SegNum].append("AllGroups")
				AllGroupsList.append(SegNum)
	if Verbose == True:
		print("In total, %d segments were among the best in at least one of the six categories.\n" % (len(SegScoreDict)))
		print("The following segments had all groups present: %s\n" % (", ".join([str(SegNum) for SegNum in AllGroupsList])))
	#now figuring out which of the segments were among the best in the greatest number of categories
	ChooseSegDict = defaultdict(list)#ChooseSegDict[NumCats] = list of SegNums that were the best in that many categories
	for SegNum in sorted(SegScoreDict.keys()):#just sorting this to make the output nicer
		if Verbose == True: print("Segment %d was among the best in the categories %s.\n" % (SegNum, ", ".join(SegScoreDict[SegNum])))
		ChooseSegDict[len(SegScoreDict[SegNum])].append(SegNum)
	ChooseSegSorted = sorted(ChooseSegDict.keys())#the sorted list of NumCats(ranging from 1 to 5, because there are at most 5 categories in which a segment can be among the best)
	BestSegList = ChooseSegDict[ChooseSegSorted[-1]]#the list of the segments that were among the best in the greatest number of categories
	#If one segment is the best in more categories, choose it.
	if len(BestSegList) == 1:
		BestSeg = BestSegList[0]
		if Verbose == True: print("The best segment is clearly %d, as it is the only segment that is among the best in %d or more categories!\n" % (BestSeg, ChooseSegSorted[-1]))
	#If not, we need to look more closely.
	else:
		#see which other segments are close for total overlap length:
		BestOVList = [ ]
		for SegNum in sorted(SegDictIn.keys()):
			if SegDictIn[SegNum]['TotOVLen'] >= GSTotOV:
				BestOVList.append(SegNum)
		if len(BestOVList) == 1:
			if Verbose == True: print("Segment %d is preferred, as it had the greatest mean overlap length.\n" % (BestOVList[0]))
			BestSeg = BestOVList[0]
		#but if there are more than one of these, looking at the rest of the categories:
		else:
			#First, try to find the absolute best segment, to see if there is a winner.
			#(I could also do something to rate the various categories.)
			#repeating the above process, but using the absolute highest values for everything, instead of the highest minus some interval
			BestScoreDict = defaultdict(list)#equivalent of SegScoreDict
			for SegNum in sorted(BestOVList):
				if SegDictIn[SegNum]['NumContigs'] == MostContigs:
					BestScoreDict[SegNum].append("NumContigs")
				if SegDictIn[SegNum]['NumInds'] == MostInds:
					BestScoreDict[SegNum].append("NumInds")
				if SegDictIn[SegNum]['MeanOVLength'] == MostOV:
					BestScoreDict[SegNum].append("MeanOVLength")
				if SegDictIn[SegNum]['MedCLength'] == LongestC:
					BestScoreDict[SegNum].append("MedCLength")
				if SegDictIn[SegNum]['TotOVLen'] == MostTotOV:
					BestScoreDict[SegNum].append("TotOVLen")
				if SegDictIn[SegNum]['NumGroups'] == TotalGroups:
					BestScoreDict[SegNum].append("AllGroups")
			if Verbose == True: print("In total, %d of these segments segments were the best in at least one of the six categories.\n" % (len(BestScoreDict)))
			ChooseBestSegDict = defaultdict(list)#This is the equivalent of ChooseSegDict
			for SegNum in sorted(BestScoreDict.keys()):
				if Verbose == True: print("Segment %d was the best overall in the categories %s.\n" % (SegNum, ", ".join(BestScoreDict[SegNum])))
				ChooseBestSegDict[len(BestScoreDict[SegNum])].append(SegNum)
			ChooseBestSegSorted = sorted(ChooseBestSegDict.keys())#equivalent to ChooseSegSorted
			WinningSegList = ChooseBestSegDict[ChooseBestSegSorted[-1]]
			#Hopefully this process will result in one segment being chosen.
			if len(WinningSegList) == 1:
				BestSeg = WinningSegList[0]
				if Verbose == True: print("The best segment is clearly %d, as it is the only segment that is the absolute best in %d or more categories!\n" % (BestSeg, ChooseBestSegSorted[-1]))
			#But, if not, the segment with the longest total overlapping length will be chosen.  Or maybe do that to start out with??
			#It would be possible to instead rate categories so that the sequence with the most individuals would be chosen, and if they were both equal, go to another category, etc.
			else:
				BestSeg = random.choice(WinningSegList)
				if Verbose == True:
					print("There were %d segments that were the absolute best in %d categories: %s.\n" % (len(WinningSegList), ChooseBestSegSorted[-1], ", ".join([str(Num) for Num in WinningSegList])))
					print("Of these segments, segment %d was randomly chosen to be used further.\n" % (BestSeg))
	if BestSeg == 0:
		sys.exit("ERROR!!!  A segment could not be selected for further analysis.")
	return BestSeg
	#This is SelectedSeg

#SegAlMaker makes a new alignment that only contains the sequences that have significant overlap with a given segment.
#the GContigList can either be "all", if all contigs for that segment should be added, or a list, if that is not true.
def SegAlMaker(SegNum, SegDictIn, GContigList, OldAlign):
	#go through the old alignment
	if GContigList == "all":
		NewAlign = AlignmentSubsetting(OldAlign,SegDictIn[SegNum]['ContigList'])
	else:
		ListTemp = list(set(SegDictIn[SegNum]['ContigList']).intersection(GContigList))
		NewAlign = AlignmentSubsetting(OldAlign,ListTemp)
		if NewAlign == "":
			print("No contigs were found for segment %d.\n" % (SegNum))
	if Verbose == True: print("%d contigs were written to the new alignment for segment %d.\n" % (len(NewAlign), SegNum))
	return NewAlign
	#This is SegmentAlign or SegmentGroupAlign

#HeaderDDPrinting
#This prints a dictionary into a file with a header.  The dictionary should have 4 levels of keys:
#The first and second keys are the first and second columns; the third key is used to sort the fourth keys; the fourth key is in the header.
#The values are in the third and remaining columns.
def HeaderDDPrinting(DictTemp, Key1Name, Key2Name, Key3List, FileName):
	OutList = [Key1Name+"\t"+Key2Name+"\t"+"\t".join(Key3List)+"\n"]
	for Key1 in sorted(DictTemp.keys()):
		for Key2 in sorted(DictTemp[Key1].keys()):
			Line = Key1+"\t"+Key2
			for Key3 in Key3List:
				try:
					Line += "\t"+str(DictTemp[Key1][Key2][Key3])
				except KeyError:
					Line += "\t"
			OutList.append(Line+"\n")
	OutFileWriting(FileName, OutList)
	
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

#########################################################################################################

#Figure out which group has which paralog for which locus
ParalogDict = DListDictFromFile(PFileName, 0, 1, 2, 3, ",")
IndDict = DictFromFile(IndDictFileName, 1, 0)
ParalogGroupDict = defaultdict(dict)
#and making a new Group entry for "all", that contains all of the contigs for that paralog:
for Locus in ParalogDict:
	for Paralog in ParalogDict[Locus]:
		ParalogGroupDict[Locus][Paralog] = defaultdict(dict)
		ParalogDict[Locus][Paralog]['all'] = [ ]
		for Group in ParalogDict[Locus][Paralog]:
			if Group != "all":
				ListTemp = ParalogDict[Locus][Paralog][Group]
				ParalogDict[Locus][Paralog]['all'] += ListTemp
				for ContigName in ListTemp:
					ParalogGroupDict[Locus][Paralog][ContigName] = Group

#make the dictionary of the taxonomic groups and their shortened names
#GroupDict[Group] = GroupAbb
GroupDict = DictFromFile(GFileName, 0, 1)
for Group in GroupDict:
	if GroupDict[Group] == "none":
		GroupDict[Group] = ""



ParalogInfoList = ["Paralog\tSegment_Number\tNumber_of_Contigs\tSegment_Length\tList_of_Groups\n"]
ISIDict = defaultdict(dict)#ISIDict[Locus][Paralog][Ind]=NumContigs
IGDict = defaultdict(list)
for Locus in ParalogDict:
	for Paralog in ParalogDict[Locus]:
		ISIDict[Locus][Paralog] = defaultdict(int)
		InFileName = InFolder+InFilePre+Paralog+"_allseqs_al.fa"
		#reading the alignment with all sequences (both our contigs and the backbone sequences)
		PAlignAllSeqs = AlignIO.read(open(InFileName), "fasta")
		PAlignGaps = AlignmentSubsetting(PAlignAllSeqs, ParalogDict[Locus][Paralog]['all'])
		if Verbose == True: print("The alignment for paralog %s has %d sequences, is %d nucleotides long with gaps, and includes the following groups: %s.\n" % (Paralog, len(PAlignGaps), PAlignGaps.get_alignment_length(), ", ".join([Group for Group in ParalogDict[Locus][Paralog].keys() if Group not in ['all']])))
		#removing the gaps from the alignment
		PAlign = AlignmentGapRemoving(PAlignGaps)
		OutFileName = InFolder+InFilePre+Paralog+"align_temp_onlycontigs.fa"
		AlignIO.write(PAlign, OutFileName, "fasta")
		if Verbose == True: print("With no gaps, the alignment has %d sequences and is %d nucleotides long.\n" % (len(PAlign), PAlign.get_alignment_length()))
		#figuring out where each contig is in the alignment
		(CInfoDict, IndSeqDict, OverDict) = AlInfoFinder(PAlign, "-", 0.8)
		#dividing the alignment into segments by looking for places where a lot of contigs start or end
		AlSegDict = BreakFinder(CInfoDict, 0.02)
		#figuring out which contigs belong to each segment
		(SegMemDict, IndSegDict) = SegMembershipFinder(AlSegDict, CInfoDict, 0.8, ParalogGroupDict[Locus][Paralog])
		#selecting the best segment based on number of contigs, number of individuals, length of sequences in that segment, and length of contigs that include that segment
		SelectedSeg = SegSelector(SegMemDict, len(ParalogDict[Locus][Paralog].keys())-1)
		#I want two versions of the alignment.  One with each group separately and one with everything together.
		#first for everything together:
		#making a dictionary that includes only the contigs that belong to that segment
		SegmentAlign = SegAlMaker(SelectedSeg, SegMemDict, "all", PAlign)
		#adding those individuals to the ISIDict
		for record in SegmentAlign:
			Ind = record.id.split("-")[0]
			ISIDict[Locus][Paralog][Ind] += 1
			Group = IndDict[Ind]
			IGDict[Group].append(Ind)
		#writing the fasta file of those contigs
		OutFileName = OutFolder+OutFilePre+Paralog+".fa"
		AlignIO.write(SegmentAlign, OutFileName, "fasta")
		print("The alignment of %d sequences of segment %d of paralog %s was written to the file %s.\n" % (len(SegmentAlign), SelectedSeg, Paralog, OutFileName))
		sys.stderr.write("The alignment of %d sequences of segment %d of paralog %s was written to the file %s.\n" % (len(SegmentAlign), SelectedSeg, Paralog, OutFileName))
		#then for each group separately:
		for Group in ParalogDict[Locus][Paralog]:
			if Group != "all":
				SegmentAlignGroup = SegAlMaker(SelectedSeg, SegMemDict, ParalogDict[Locus][Paralog][Group], PAlign)
				OutFileName = SeqFilePath+Group+"/"+SeqFolderPre+GroupDict[Group]+"/"+SeqFolderPost+SeqFilePre+Paralog+"seg_"+str(SelectedSeg)+".fa"
				if (SeqFolderPre == "") and (GroupDict[Group] == ""):
					OutFileName = SeqFilePath+Group+"/"+SeqFolderPost+SeqFilePre+Paralog+"seg_"+str(SelectedSeg)+".fa"
				AlignIO.write(SegmentAlign, OutFileName, "fasta")
		print("Individual alignments for the %d groups were written to files with names such as %s.\n" % (len(ParalogDict[Locus][Paralog].keys())-1, OutFileName))
		sys.stderr.write("Individual alignments for the %d groups were written to files with names such as %s.\n" % (len(ParalogDict[Locus][Paralog].keys())-1, OutFileName))
	#adding information about this paralog to the list of paralogs
	Line = Paralog+"\t"+str(SelectedSeg)+"\t"+str(len(SegmentAlign))+"\t"+str(SegMemDict[SelectedSeg]['SegEnd']-SegMemDict[SelectedSeg]['SegStart']+1)+"\t"+", ".join([Item for Item in ParalogDict[Locus][Paralog].keys() if Item != "all"])+"\n"
	ParalogInfoList.append(Line)

OutFileName = OutFolder+OutFilePre+"Paralog_Info.txt"
OutFileWriting(OutFileName, ParalogInfoList)

for Group in IGDict:
	ListTemp = list(set(IGDict[Group]))
	IGDict[Group] = ListTemp
IndList = [ ]
for Group in sorted(IGDict.keys()):
	for Ind in sorted(IGDict[Group]):
		IndList.append(Ind)
HeaderDDPrinting(ISIDict, "Locus", "Paralog", IndList, OutFolder+OutFilePre+"Seqs_per_Paralog.txt")
