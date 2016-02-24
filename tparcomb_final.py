#! /usr/bin/env python

#tparcomb_final.py version 1.0 31 Dec. 2015 Abby Moore
#This script takes the results of tparcomb_combiner.py and deals with the loci that had
#multiple sequences for the same paralog for some individuals.
#Once it has fixed these loci, it makes final versions of two information files.
#version 1.1 23 Feb. 2016
#modified so that it aligns the sequences and makes trees in order from largest to smallest sequence file
#**************This has not been tested!!********************

#format of RedoList.txt (has a header line, starting with "Locus" that provides no useful information)
'''
Locus: [0]
Paralog: [1]
Group: [2]
list of Ind Names: [3]
'''
#format of Ind_Seq_Info.txt (has a header line, starting with "Locus" that identifies the individuals)
'''
Locus: amk [0]
Paralog: Ambig_amk_none [1]
individual information: 12 seqs, overlap: 630, 59 ambigs, length: 657, using: no [2 and all subsequent columns]
Columns are blank if there are no sequences for that paralog from that group
'''
#format of Locus_Paralog_List.txt (no header)
'''
Locus: amk [0]
Paralog: amk1 [1]
'''

import sys
from collections import defaultdict
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

Usage = '''
tparcomb_final.py combines multiple sequences from the same paralog for some individuals
that were uncovered by tparcomb_combiner.py and then writes final versions of the 
alignments and other files.
[input folder]
[prefix for input files]
[output folder]
[prefix for output files]
[path to folder where scripts are found, or "none", if none]
[the folder for the backbone alignments and trees]
[the prefix for the backbone alignments and trees]
[the suffix for the backbone alignments]
[outgroup dictionary]
[dictionary saying which group each individual belongs to]
[number of cores for parallelization]
'''

'''
tparcomb_final.py InFilePath InFileGroupList IndDictFileName OutGroupDictFN ContigFolder ContigFilePre OutFolder OutFilePre AlFolder AlFilePre AlFilePost ScriptPath OutGroupDictFN IndDictFileName NCores
/gpfs/scratch/ajm3/eedwards/scripts/tparcomb_final.py /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_combined Ln1tcb_ /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_final Ln1tfi_ /gpfs/scratch/ajm3/eedwards/scripts/ /gpfs/scratch/ajm3/eedwards/general/combined_trees/ new_ _al /gpfs/scratch/ajm3/eedwards/general/outgroup_list_new.txt /gpfs/scratch/ajm3/eedwards/general/Ln1_inds2.txt 4
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
OutFolder = sys.argv[3]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[4]
if OutFilePre == "none":
	OutFilePre = ""
ScriptPath = sys.argv[5]
if ScriptPath[-1] != "/":
	ScriptPath += "/"
if ScriptPath == "none/":
	ScriptPath = ""
AlFolder = sys.argv[6]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[7]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[8]
if AlFilePost == "none":
	AlFilePost = ""
OutGroupDictFN = sys.argv[9]
IndDictFileName = sys.argv[10]
NCores = sys.argv[11]

#SeqFileReading reads a sequence file and puts the sequences in a dictionary.
def SeqFileReading(FileName, SeqFormat):
	DictTemp = { }#DictTemp[SeqName] = Seq
	InFile = open(FileName, 'rU')
	for record in SeqIO.parse(InFile, SeqFormat):
		DictTemp[record.id] = str(record.seq)
	InFile.close()
	return DictTemp
	#This is various sequence dictionaries

#SeqFileWriting writes sequence files from dictionaries.
def SeqFileWriting(FileName, SDict, SeqFormat):
	OutFile = open(FileName, 'w')
	for SeqName in SDict:
		Record1 = SeqRecord(seq=Seq(SDict[SeqName], IUPAC), id = SeqName, description = "")
		SeqIO.write(Record1, OutFile, SeqFormat)
	OutFile.close()

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
	#This is OutGroupDict
	
#HeaderDDMaking
#Key1Name is the name of the first word of the header line, which is in the column that will be the first key when reading the other lines.
#GroupDict is the dictionary used to classify the groups.
#This is modified from tparcomb_combiner.py, where it was not a function
def HeaderDDMaking(FileName, Key1Name):
	DictTemp = defaultdict(dict)
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		if Line[0] == Key1Name:
			#the line is the list of keys for which individuals are in which positions
			KeyList = Line
			KeyRange = range(2,len(KeyList))
		else:
			Key1 = Line[0]
			Key2 = Line[1]
			if Line[2:] == [""]*(len(Line[2:])):
				print("No information for %s: %s.\n" % (Key1, Key2))
			else:
				for KeyPos in KeyRange:
					Key3= KeyList[KeyPos]
					MyValue = Line[KeyPos]
					try:
						DictTemp[Key1][Key2][Key3] = MyValue
					except KeyError:
						DictTemp[Key1][Key2] = defaultdict(dict)
						DictTemp[Key1][Key2][Key3] = MyValue
	InFile.close()
	return (DictTemp, KeyList[2:])
	#These are ISIDict and IndList.

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
			if (Line == Key1+"\t"+Key2+("\t"*len(Key3List))) or (Line == Key1+"\t"+Key2+("0\t"*len(Key3List))):
				print("The line for %s: %s is now empty.\n" % (Key1, Key2))
			else:
				OutList.append(Line+"\n")
	OutFileWriting(FileName, OutList)


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

#########################################################################################################################

#First, read the RedoList file.
RedoList = defaultdict(dict)#RedoList[Locus][Paralog][Group] = list of individuals
InFileName = InFolder+InFilePre+"RedoList.txt"
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	#if this is not the header line:
	if Line[0] != "Locus":
		Locus = Line[0]
		Paralog = Line[1]
		Group = Line[2]
		IndList = Line[3].split(', ')
		try:
			RedoList[Locus][Paralog][Group] = IndList
		except KeyError:
			RedoList[Locus][Paralog] = defaultdict(list)
			RedoList[Locus][Paralog][Group] = IndList
InFile.close()

#making the out group dictionary
OutGroupDict = DictFromFile(OutGroupDictFN, 0, 1)
IndDict = DictFromFile(IndDictFileName, 1, 0)

#The question is whether I should read all of the other files into dictionaries, so I can change them as I go, or whether
#I should deal with the sequences first and then change the files as I read them in.

NewSeqDict = defaultdict(dict)#NewSeqDict[Locus][Paralog][Ind]['SeqName']/['Seq']/['Fate']
RedoneDict = defaultdict(list) #RedoneDict[Locus] = list of ParalogNames
for Locus in RedoList:
	NewLocusSeqs = { }
	#I will be using _allseqs_al.fa as a source for the sequences
	LocusSeqDict = SeqFileReading(InFolder+InFilePre+Locus+"_allseqs_al.fa", "fasta")
	#But other sequence files will need to be modified as well.
	BestSeqDict = SeqFileReading(InFolder+InFilePre+Locus+"_allbest_al.fa", "fasta")
	AllCombSeqDict = SeqFileReading(InFolder+InFilePre+Locus+"_combined_all.fa", "fasta")
	BestCombSeqDict = SeqFileReading(InFolder+InFilePre+Locus+"_combined_best.fa", "fasta")
	for Paralog in RedoList[Locus]:
		NewSeqDict[Locus][Paralog] = defaultdict(dict)
		ParalogName = Paralog.split("_P")[0]
		ParalogSeqDict = SeqFileReading(InFolder+InFilePre+Locus+"_"+ParalogName+".fa", "fasta")
		for Group in RedoList[Locus][Paralog]:#I am not sure that we care about the group.
			for Ind in RedoList[Locus][Paralog][Group]:
				IndSeqList = [ ]
				NameTemp = Ind+"."+Paralog
				for SeqName in LocusSeqDict:
					if SeqName[:len(NameTemp)] == NameTemp:
						IndSeqList.append(SeqName)
				(ConSeq, NumAmbig, AmbigSitesList, NumSeqs, Overlap) = ConSeqMaker(LocusSeqDict,IndSeqList)
				#keep the sequence, if there are few enough ambiguities
				if NumAmbig < 5:
					#First, add everything about that sequence to the dictionary
					NewSeqDict[Locus][Paralog][Ind]['Fate'] = 'combine'
					NewSeqDict[Locus][Paralog][Ind]['Seq'] = ConSeq
					NumSeqs = 0
					NumAmbigs = NumAmbig
					for SeqName in IndSeqList:
						SeqNameSplit = SeqName.split(".")
						NumSeqs += int(SeqNameSplit[2][:-4])
						NumAmbigs += int(SeqNameSplit[3][:-5])
					NewSeqName = ".".join([Ind, Paralog, str(NumSeqs)+"seqs", str(NumAmbigs)+"ambig", "len"+str(len(ConSeq))])
					NewSeqDict[Locus][Paralog][Ind]['SeqName'] = NewSeqName
					NewSeqDict[Locus][Paralog][Ind]['SeqLen'] = len(ConSeq)
					NewSeqDict[Locus][Paralog][Ind]['NumAmbigs'] = NumAmbigs
					NewSeqDict[Locus][Paralog][Ind]['NumSeqs'] = NumSeqs
					NewSeqDict[Locus][Paralog][Ind]['Overlap'] = Overlap
					NewSeqDict[Locus][Paralog][Ind]['OldSeqs'] = IndSeqList
					#Then, remove both of the old sequences from the alignments
					for SeqName in IndSeqList:
						del LocusSeqDict[SeqName]
						del ParalogSeqDict[SeqName]
						del AllCombSeqDict[SeqName]
						try:
							del BestSeqDict[SeqName]
							del BestCombSeqDict[SeqName]
						except KeyError:
							"do nothing"
					#And add the new sequences to the dictionaries of unaligned sequences
					AllCombSeqDict[NewSeqName] = ConSeq
					BestCombSeqDict[NewSeqName] = ConSeq
					NewLocusSeqs[NewSeqName] = ConSeq
					#and the paralog dictionary
					ParalogSeqDict[NewSeqName] = ConSeq
					RedoneDict[Locus].append(ParalogName)
				else:
					"do nothing"
					NewSeqDict[Locus][Paralog][Ind]['Fate'] = 'separate'
				#If there are too many ambiguities, then I guess we just keep the sequences separate, because they are probably different.
		#write the new paralog file (the old file with the duplicate sequences removed and the new paralog sequences added)
		SeqFileWriting(OutFolder+OutFilePre+Locus+"_"+ParalogName+".fa", ParalogSeqDict, "fasta")
	#write the new locus files
	SeqFileWriting(OutFolder+OutFilePre+Locus+"_allseqs.fa", LocusSeqDict, "fasta")
	SeqFileWriting(OutFolder+OutFilePre+Locus+"_allbest.fa", BestSeqDict, "fasta")
	SeqFileWriting(OutFolder+OutFilePre+Locus+"_combined_all.fa", AllCombSeqDict, "fasta")
	SeqFileWriting(InFolder+InFilePre+Locus+"_combined_best.fa", BestCombSeqDict, "fasta")
	#write the file of sequences to be aligned to the new locus files
	SeqFileWriting(OutFolder+OutFilePre+Locus+"_newseqs.fa", NewLocusSeqs, "fasta")
	print("The updated sequence files for locus %s were written.\n" % (Locus))
	sys.stderr.write("The updated sequence files for locus %s were written.\n" % (Locus))

#############################################################################
#Updating the various information files
#Ambigs_per_Seq.txt is not going to be easy to figure out at this point.  Maybe we don't need a final version of it??
#Contig_Fates.txt really doesn't change much, so that can just be left as is.
#the Seq_Fates files aren't that interesting, so we don't really need them.
#Seqs_per_Locus.txt is also not that interesting, because it's grouped by group, not by individual.
#perc_ambig.txt needs to be updated but won't be hard.  I am also not sure that this is very interesting, given that I have all of those different
#sub-paralogs.
#Locus_Paralog_List.txt doesn't change, so doesn't need to be updated.

#Ind_Seq_Info.txt
#reading in the old file:
#(ISIDict, IndList) = HeaderDDMaking(InFolder+InFilePre+"Ind_Seq_Info.txt", "Locus", IndDict)
(ISIDict, IndList) = HeaderDDMaking(InFolder+InFilePre+"Ind_Seq_Info.txt", "Locus")
#getting information on the number of sequences per paralog from the ISIDict:
IPDict = defaultdict(dict)#IPDict[Locus][Paralog][Ind] = NumParalogs
#without Groups
for Locus in ISIDict:
	for Paralog in ISIDict[Locus]:
		ParalogTemp = Paralog.split("_P")[0]
		for Ind in ISIDict[Locus][Paralog]:
			#if Ind == "Mollugo_verticillata_56":
			#	print("%s: %s: %s" % (Locus,Paralog, ISIDict[Locus][Paralog][Ind]))
			#if I am using the sequence directly (and not the contig), then add it to dictionary of information about good sequences
			if ISIDict[Locus][Paralog][Ind][-3:] == "yes":
				try:
					IPDict[Locus][ParalogTemp][Ind] += 1
				except KeyError:
					IPDict[Locus][ParalogTemp] = defaultdict(int)
					IPDict[Locus][ParalogTemp][Ind] = 1
#substituting the new information:
for Locus in NewSeqDict:
	for Paralog in NewSeqDict[Locus]:
		for Ind in NewSeqDict[Locus][Paralog]:
			if NewSeqDict[Locus][Paralog][Ind]['Fate'] == 'combine':
				NewLine = "round comb, "+str(NewSeqDict[Locus][Paralog][Ind]['NumSeqs'])+" seqs, overlap: "+str(NewSeqDict[Locus][Paralog][Ind]['Overlap'])+", "+str(NewSeqDict[Locus][Paralog][Ind]['NumAmbigs'])+" ambigs, length: "+str(NewSeqDict[Locus][Paralog][Ind]['SeqLen'])+", using: yes"
				ISIDict[Locus][Paralog][Ind] = NewLine
			elif NewSeqDict[Locus][Paralog][Ind]['Fate'] == 'separate':
				#at this point, I think I will just leave that in the ISIDict.
				#It is good information that they may or may not actually be separate paralogs, anyway.
				ParalogTemp = Paralog.split("_P")[0]
				IPDict[Locus][ParalogTemp][Ind] += 1


HeaderDDPrinting(ISIDict, "Locus", "Paralog", IndList, OutFolder+OutFilePre+"Ind_Seq_Info.txt")
HeaderDDPrinting(IPDict, "Locus", "Paralog", IndList, OutFolder+OutFilePre+"Seqs_per_Paralog.txt")

#making new sequence files with all sequences that are over 100bp
OutLocusDict = defaultdict(list)
for Locus in ISIDict:
	NumSeqs = 0
	if Locus in NewSeqDict.keys():
		InFileName = OutFolder+OutFilePre+Locus+"_allseqs.fa"
	else:
		InFileName = InFolder+InFilePre+Locus+"_allseqs_al.fa"
	OutFileName = OutFolder+OutFilePre+Locus+"_allover100.fa"
	InFile = open(InFileName, 'rU')
	OutFile = open(OutFileName, 'w')
	for record in SeqIO.parse(InFile, "fasta"):
		NumSeqs += 1
		if len(str(record.seq).replace("-","")) > 100:
			SeqIO.write(record, OutFile, "fasta")
	InFile.close()
	OutFile.close()
	OutLocusDict[NumSeqs].append(Locus)
print("%d sequence files that only contain the sequences that are longer than 100 bp were written, with names such as %s.\n" % (len(ISIDict.keys()), OutFileName))
sys.stderr.write("%d sequence files that only contain the sequences that are longer than 100 bp were written, with names such as %s.\n" % (len(ISIDict.keys()), OutFileName))


#updating pardict
#reading in the old pardict:
InFileName = InFolder+InFilePre+"pardict_new.txt"
InFile = open(InFileName, 'rU')
PDDict = defaultdict(dict)#PDDict[Locus][Paralog][SeqName] = Genus
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	Locus = Line[0]
	SeqName = Line[1]
	Paralog = Line[2]
	Genus = Line[3]
	try:
		PDDict[Locus][Paralog][SeqName] = Genus
	except KeyError:
		PDDict[Locus][Paralog] = defaultdict(dict)
		PDDict[Locus][Paralog][SeqName] = Genus
InFile.close()
#adding the changed sequences to the pardict:
for Locus in NewSeqDict:
	for Paralog in NewSeqDict[Locus]:
		ParalogTemp = Paralog.split("_P")[0]
		for Ind in NewSeqDict[Locus][Paralog]:
			if NewSeqDict[Locus][Paralog][Ind]['Fate'] == 'combine':
				for OldSeq in NewSeqDict[Locus][Paralog][Ind]['OldSeqs']:
					del PDDict[Locus][ParalogTemp][OldSeq]
				SeqName = NewSeqDict[Locus][Paralog][Ind]['SeqName']
				Genus = Ind.split("_")[0]
				PDDict[Locus][ParalogTemp][SeqName] = Genus
#writing the new pardict:
OutList = [ ]
for Locus in sorted(PDDict.keys()):
	for Paralog in sorted(PDDict[Locus].keys()):
		for SeqName in sorted(PDDict[Locus][Paralog].keys()):
			Line = '\t'.join([Locus, SeqName, Paralog, PDDict[Locus][Paralog][SeqName]])
			OutList.append(Line+"\n")
OutFileName = OutFolder+OutFilePre+'pardict.txt'
OutFileWriting(OutFileName, OutList)

#############################################################################
#Writing a script to transfer the other files, align the new sequences to the old files, and make trees.

InFileName = InFolder+InFilePre+"Locus_Paralog_List.txt"
InFile = open(InFileName, "rU")
LPDict = defaultdict(list)
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	Locus = Line[0]
	Paralog = Line[1]
	if len(Paralog.split('_P')) == 1:
		LPDict[Locus].append(Paralog)
InFile.close()

OutLocusList = [ ]
for NumSeqs in sorted(OutLocusDict.keys(), reverse=True):
	if NumSeqs != 0:
		OutLocusList += OutLocusDict[NumSeqs]

OutList1 = ['#! /bin/bash\n']
OutFileName1 = OutFolder+OutFilePre+"Locus_Analysis_Subscript1.sh"
OutList2 = ['#! /bin/bash\n']
OutFileName2 = OutFolder+OutFilePre+"Locus_Analysis_Subscript2.sh"
OutList3 = ['#! /bin/bash\n']
OutFileName3 = OutFolder+OutFilePre+"Locus_Analysis_Script.sh"
for Locus in OutLocusList:
	#if we don't need to add any new sequences to that locus
	if (Locus in RedoneDict.keys()) == False:
		Line = 'cp '+InFolder+InFilePre+Locus+"_allbest_al.fa "+OutFolder+OutFilePre+Locus+"_allbest_al.fa "
		Line += "&& "+ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allbest_al.fa\n"
		Line += 'cp '+InFolder+InFilePre+Locus+"_allseqs_al.fa "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa "
		Line += "&& "+ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		Line += "mv "+OutFolder+OutFilePre+Locus+"_allover100.fa "+OutFolder+OutFilePre+Locus+"_allover100_al.fa  "
		Line += "&& "+ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allover100_al.fa\n"
		OutList1.append(Line)
		#and the paralogs
		for Paralog in LPDict[Locus]:
			if Paralog[:5] != "Ambig":
				Line = 'cp '+InFolder+InFilePre+Locus+"_"+Paralog+".fa "+OutFolder+OutFilePre+Locus+"_"+Paralog+".fa "
				Line += "&& mafft --addfragments "+OutFolder+OutFilePre+Locus+"_"+Paralog+".fa --quiet --thread -1 "+AlFolder+AlFilePre+Paralog+AlFilePost+".fa > "+OutFolder+OutFilePre+Locus+"_"+Paralog+"_al.fa\n"
				OutList1.append(Line)
	#if we do need to add new sequences to that locus
	else:
		Line = "mafft --addfragments "+OutFolder+OutFilePre+Locus+"_newseqs.fa --quiet --thread -1 "+OutFolder+OutFilePre+Locus+"_allseqs.fa > "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa "
		Line += "&& "+ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa\n"
		Line += "mafft --addfragments "+OutFolder+OutFilePre+Locus+"_newseqs.fa --quiet --thread -1 "+OutFolder+OutFilePre+Locus+"_allbest.fa > "+OutFolder+OutFilePre+Locus+"_allbest_al.fa "
		Line += "&& "+ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allbest_al.fa\n"
		Line += "mafft --addfragments "+OutFolder+OutFilePre+Locus+"_newseqs.fa --quiet --thread -1 "+OutFolder+OutFilePre+Locus+"_allover100.fa > "+OutFolder+OutFilePre+Locus+"_allover100_al.fa "
		Line += "&& "+ScriptPath+"fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+"_allover100_al.fa\n"
		OutList1.append(Line)
		for Paralog in LPDict[Locus]:
			if Paralog[:5] != "Ambig":
				if (Paralog in RedoneDict[Locus]) == False:
					Line = 'cp '+InFolder+InFilePre+Locus+"_"+Paralog+".fa "+OutFolder+OutFilePre+Locus+"_"+Paralog+".fa "
					Line += "&& mafft --addfragments "+OutFolder+OutFilePre+Locus+"_"+Paralog+".fa --quiet --thread -1 "+AlFolder+AlFilePre+Paralog+AlFilePost+".fa > "+OutFolder+OutFilePre+Locus+"_"+Paralog+"_al.fa\n"
					OutList1.append(Line)
				else:
					Line = "mafft --addfragments "+OutFolder+OutFilePre+Locus+"_"+Paralog+".fa --quiet --thread -1 "+AlFolder+AlFilePre+Paralog+AlFilePost+".fa > "+OutFolder+OutFilePre+Locus+"_"+Paralog+"_al.fa\n"
					OutList1.append(Line)
	#These trees are just going to be used to classify the various paralogs, so I do not need to have bootstrap values
	Line = "raxmlHPC  -f d -s "+OutFolder+OutFilePre+Locus+"_allbest_al.phy -n "+OutFilePre+"ab_"+Locus+" -m GTRCAT -p 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
	#Line += "raxmlHPC -f d -s "+OutFolder+OutFilePre+Locus+"_allseqs_al.phy -n "+OutFilePre+Locus+" -m GTRCAT -p 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
	Line += "raxmlHPC -f d -s "+OutFolder+OutFilePre+Locus+"_allover100_al.phy -n "+OutFilePre+"o100_"+Locus+" -m GTRCAT -p 1234 -o "+OutGroupDict[Locus]+" -w "+OutFolder+"\n"
	OutList2.append(Line)
OutFileWriting(OutFileName1, OutList1)
Line = "chmod u+x "+OutFileName1+"\n"
Line += "cat "+OutFileName1+" | parallel --jobs "+NCores+" --joblog "+OutFolder+OutFilePre+"parallel_log1.log\n"
OutList3.append(Line)
OutFileWriting(OutFileName2, OutList2)
Line = "chmod u+x "+OutFileName2+"\n"
Line += "cat "+OutFileName2+" | parallel --jobs "+NCores+" --joblog "+OutFolder+OutFilePre+"parallel_log2.log\n"
OutList3.append(Line)
OutFileWriting(OutFileName3, OutList3)
print("The shell script for the sequences was written to %s.\n" % (OutFileName3))
sys.stderr.write("The shell script for the sequences was written to %s.\n" % (OutFileName3))
