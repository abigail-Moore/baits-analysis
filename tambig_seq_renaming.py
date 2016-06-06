#! /usr/bin/env python

#tambig_seq_renaming.py version 1.0 7 Jan. 2016 Abby Moore
#This script renames contigs according to an input file of the following form:
#Ambig_Paralog_List.txt
'''
Locus [0]: asp
OldParalog [1]: Ambig_asp_3
GroupList [2]: All or a comma-delimited list of groups
NewParalog [3]: asp31
'''

import sys
from collections import defaultdict
import os #We need to be able to talk to the operating system
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

Usage = '''
tambig_seq_renaming.py takes ambiguous sequences that were not classified 
according to paralog and adds them to the user-chosen paralogs.
[input folder]
[prefix for input files]
[contig folder]
[prefix for contig files]
[combined sequence folder]
[prefix for combined sequences]
[output folder, or "same", if the same as the input folder]
[prefix for output files, or "same", if the same as the input file prefix]
[path to folder where scripts are found, or "none", if none]
[outgroup dictionary]
[dictionary saying which group each individual belongs to]
[number of cores for parallelization]
[file saying which paralogs need to be in which groups]
'''

'''
tambig_seq_renaming.py InFolder InFilePre CFolder CFilePre CBFolder CBFilePre OutFolder OutFilePre ScriptPath OutGroupDictFN IndDictFileName NCores NPFileName
/gpfs/scratch/ajm3/eedwards/scripts/tambig_seq_renaming.py /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_final/ Ln1tfi_ /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_contigsplit/ Ln1tcs_ /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_combined/ Ln1tcb_ /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_final2/ Ln1tfi2_ /gpfs/scratch/ajm3/eedwards/scripts/ /gpfs/scratch/ajm3/eedwards/general/outgroup_list_new.txt /gpfs/scratch/ajm3/eedwards/general/Ln1_inds2.txt 8 /gpfs/scratch/ajm3/eedwards/Ln1_testing/Ln1t_combined/Ln1tcb_Ambig_Paralog_List.txt
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 14:
	sys.exit("ERROR!  This script requires 13 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
if InFolder[-1] != "/":
	InFolder += "/"
InFilePre = sys.argv[2]
if InFilePre == "none":
	InFilePre = ""
CFolder = sys.argv[3]
if CFolder[-1] != "/":
	CFolder += "/"
CFilePre = sys.argv[4]
if CFilePre == "none":
	CFilePre = ""
CBFolder = sys.argv[5]
if CBFolder[-1] != "/":
	CBFolder += "/"
CBFilePre = sys.argv[6]
if CBFilePre == "none":
	CBFilePre = ""
OutFolder = sys.argv[7]
if OutFolder == "same":
	OutFolder = InFolder
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[8]
if OutFilePre == "none":
	OutFilePre = ""
ScriptPath = sys.argv[9]
if ScriptPath[-1] != "/":
	ScriptPath += "/"
if ScriptPath == "none/":
	ScriptPath = ""
OutGroupDictFN = sys.argv[10]
IndDictFileName = sys.argv[11]
NCores = sys.argv[12]
NPFileName = sys.argv[13]

#####################################################################################################################################

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
	#This is OutGroupDict and IndDict

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
	#This is NPDict
	
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
			if (Line == Key1+"\t"+Key2+("\t"*len(Key3List))) or (Line == Key1+"\t"+Key2+("\t0"*len(Key3List))):
				print("The line for %s: %s is now empty.\n" % (Key1, Key2))
			else:
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


##############################################################################################################################

#Read the file saying which paralogs need to be renamed in which groups
#NPDict[Locus][OldPar][NewPar] = list of Groups
NPDict = DListDictFromFile(NPFileName,0,1,3,2,",")

#Read the file saying which individuals belong to which groups
#IndDict[Ind] = Group
IndDict = DictFromFile(IndDictFileName, 1, 0)
#GroupDict[Group] = list of Inds
GroupDict = defaultdict(list)
for Ind in IndDict:
	GroupDict[IndDict[Ind]].append(Ind)

#then the outgroups
#OutGroupDict[Locus] = Outgroup
OutGroupDict = DictFromFile(OutGroupDictFN, 0, 1)

#make the Ind_Seq_Info and Seqs_per_Paralog dictionaries so I can update them.
(ISIDict, ISIIndList) = HeaderDDMaking(InFolder+InFilePre+"Ind_Seq_Info.txt", "Locus")
(IPDict, IPIndList) = HeaderDDMaking(InFolder+InFilePre+"Seqs_per_Paralog.txt", "Locus")
#and the same for the contig equivalent
(CIPDict, CIPIndList) = HeaderDDMaking(CFolder+CFilePre+"Seqs_per_Paralog.txt", "Locus")
#make the pardict dicationary:
InFileName = InFolder+InFilePre+"pardict.txt"
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
if (InFolder == OutFolder) and (InFilePre == OutFilePre):
	OutLine = "mv "+InFolder+InFilePre+"Ind_Seq_Info.txt "+InFolder+InFilePre+"Ind_Seq_Info_old.txt"
	os.popen(OutLine,'r')
	OutLine = "mv "+InFolder+InFilePre+"Seqs_per_Paralog.txt "+InFolder+InFilePre+"Seqs_per_Paralog_old.txt"
	os.popen(OutLine,'r')
	OutLine = "mv "+InFolder+InFilePre+"pardict.txt "+InFolder+InFilePre+"pardict_old.txt"
OutLine = "mv "+CFolder+CFilePre+"Seqs_per_Paralog.txt "+CFolder+CFilePre+"Seqs_per_Paralog_old.txt"
os.popen(OutLine,'r')

#Things that need to be done:
#check to see if there is already a sequence from that paralog for that individual
#If there is, do I try to combine them, or do I assume they are different, because they came out in different places??
#If I don't need to combine them, then it's easy:
#rename NewPar to NewPar_old.fa
#move the sequences from the OldPar to the NewPar files (probably best to align the new sequences to the (old) NewPar instead of aligning everything)
#make a file with the old sequence name[tab]new sequence name to
#rename allbest_al.fa, allover150_al.fa, allseqs_al.fa, RAxML_bestTree..._ab_ and RAxML_bestTree..._o150
#copy at least .fa and the RAxML files to filename_old.extension, so that I can see what happened (if InFolder and OutFolder are the same)
#update the ISIDict and IPDict by removing the old paralog (if everything has been classified) and adding entries to the new paralog.

RenamedLocusList = [ ]
ParFateDict = defaultdict(dict)
NewParDict = defaultdict(list)
for Locus in NPDict:
	LocusRenamingList = [ ]
	for OldPar in NPDict[Locus]:
		try:
			OldParSeqDict = SeqFileReading(InFolder+InFilePre+Locus+"_"+OldPar+".fa", "fasta")
		except IOError:
			OldParSeqDict = SeqFileReading(CBFolder+CBFilePre+Locus+"_"+OldPar+".fa", "fasta")
		if (InFolder == OutFolder) and (InFilePre == OutFilePre):
			OutLine = "mv "+InFolder+InFilePre+Locus+"_"+OldPar+".fa "+InFolder+InFilePre+Locus+"_"+OldPar+"_old.fa"
			os.popen(OutLine,'r')
		#Sequence names of two possible forms:
		#Anacampseros_recurvata_88.Ambig_asp_3.5seqs.1ambig.len1173
		#Portulaca_cryptopetala_49-final-NODE_1_length_1424_cov_17.3214_ID_1_exons_Ambig_asp_3
		#making a dictionary of which sequences go with which individuals in the OldPar
		OldParIndDict = defaultdict(list)
		ContigIndDict = { }
		for SeqName in OldParSeqDict:
			if len(SeqName.split("_exons")) == 2:
				#it's a contig sequence
				Ind = SeqName.split("-")[0]
				ContigIndDict[Ind] = 1
			else:
				Ind = SeqName.split(".")[0]
				ContigIndDict[Ind] = 0
			OldParIndDict[Ind].append(SeqName)
		#now, looking at each NewPar
		for NewPar in NPDict[Locus][OldPar]:
			NewParDict[Locus].append(NewPar)
			NewParSeqDict = { }
			#making a list of individuals I want to classify in the NewPar
			IndsWanted = [ ]
			for Group in NPDict[Locus][OldPar][NewPar]:
				if (Group == "All") or (Group == "all"):
					IndsWanted = OldParIndDict.keys()
				else:
					for Ind in GroupDict[Group]:
						if Ind in OldParIndDict.keys():
							IndsWanted.append(Ind)
			for Ind in IndsWanted:
				#if we are dealing with contigs:
				if ContigIndDict[Ind] == 1:
					for OldSeqName in OldParIndDict[Ind]:
						#adding the sequence to the dictionary of NewPar sequences
						NewSeqName = OldSeqName.split(OldPar)[0]+NewPar
						NewParSeqDict[NewSeqName] = OldParSeqDict[OldSeqName]
						LocusRenamingList.append(OldSeqName+"\t"+NewSeqName+"\n")
						#removing it from the dictionary of OldPar sequences
						del OldParSeqDict[OldSeqName]
						#and updating the CIPDict
						try:
							CIPDict[Locus][NewPar][Ind] = int(CIPDict[Locus][NewPar][Ind]) + 1
						except KeyError:
							CIPDict[Locus][NewPar] = defaultdict(int)
							CIPDict[Locus][NewPar][Ind] = '1'
						CIPDict[Locus][OldPar][Ind] = int(CIPDict[Locus][OldPar][Ind]) - 1
				#if we are dealing with full sequences:
				elif ContigIndDict[Ind] == 0:
					for OldSeqName in OldParIndDict[Ind]:
						#adding the sequence to the dictionary of NewPar sequences
						NewSeqName = OldSeqName.split(OldPar)[0]+NewPar+OldSeqName.split(OldPar)[1]
						NewParSeqDict[NewSeqName] = OldParSeqDict[OldSeqName]
						LocusRenamingList.append(OldSeqName+"\t"+NewSeqName+"\n")
						#removing it from the dictionary of OldPar sequences
						del OldParSeqDict[OldSeqName]
						#updating the IPDict
						try:
							IPDict[Locus][NewPar][Ind] = int(IPDict[Locus][NewPar][Ind]) + 1
						except KeyError:
							IPDict[Locus][NewPar] = defaultdict(int)
							IPDict[Locus][NewPar][Ind] = 1
						IPDict[Locus][OldPar][Ind] = int(IPDict[Locus][OldPar][Ind]) - 1
						#if there isn't an entry for the NewPar yet
						#round 3, 2 seqs, overlap: 102, 5 ambigs, length: 525, using: yes
						#Anacampseros_recurvata_88.Ambig_asp_3.5seqs.1ambig.len1173
						SplitEntry = ISIDict[Locus][OldPar][Ind].split(", ")
						NewRound = ["round "+SplitEntry[0][6:]+" mv"]
						NewRound += SplitEntry[1:]
						try:
							if ISIDict[Locus][NewPar][Ind] == "":
								ISIDict[Locus][NewPar][Ind] =  ", ".join(NewRound)
							else:
								ISIDict[Locus][NewPar][Ind] += "/"+", ".join(NewRound)
								print("Warning! %s already has a sequence for paralog %s of locus %s.\n" % (Ind, NewPar, Locus))
						except KeyError:
							ISIDict[Locus][NewPar] = defaultdict(dict)
							ISIDict[Locus][NewPar][Ind] =  ", ".join(NewRound)
						except TypeError:
							ISIDict[Locus][NewPar][Ind] =  ", ".join(NewRound)
						del PDDict[Locus][OldPar][OldSeqName]
						Genus = Ind.split("_")[0]
						try:
							PDDict[Locus][NewPar][NewSeqName] = Genus
						except KeyError:
							PDDict[Locus][NewPar] = defaultdict(dict)
							PDDict[Locus][NewPar][NewSeqName] = Genus
			#write the sequences for the new paralog
			SeqFileWriting(OutFolder+OutFilePre+Locus+"_"+NewPar+".fa", NewParSeqDict, "fasta")
		#check to see if there are any OldPar sequence left.  If so, then write a new OldPar file.
		#if all of the sequences from OldPar were classified
		if OldParSeqDict == { }:
			ParFateDict[Locus][OldPar] = "no seqs remaining"
		#if not,
		else:
			ParFateDict[Locus][OldPar] = str(len(OldParSeqDict))+" seqs remaining"
			SeqFileWriting(OutFolder+OutFilePre+Locus+"_"+OldPar+".fa", OldParSeqDict, "fasta")
	#write the renaming file for the locus
	if LocusRenamingList != [ ]:
		LocusRenamingFN = OutFolder+OutFilePre+Locus+"_renaming_key.txt"
		OutFileWriting(LocusRenamingFN, LocusRenamingList)
		RenamedLocusList.append(Locus)

#writing the information files
HeaderDDPrinting(ISIDict, "Locus", "Paralog", ISIIndList, OutFolder+OutFilePre+"Ind_Seq_Info.txt")
HeaderDDPrinting(IPDict, "Locus", "Paralog", IPIndList, OutFolder+OutFilePre+"Seqs_per_Paralog.txt")
HeaderDDPrinting(CIPDict, "Locus", "Paralog", CIPIndList, CFolder+CFilePre+"Seqs_per_Paralog.txt")
#writing the new pardict:
OutList = [ ]
for Locus in sorted(PDDict.keys()):
	for Paralog in sorted(PDDict[Locus].keys()):
		for SeqName in sorted(PDDict[Locus][Paralog].keys()):
			Line = '\t'.join([Locus, SeqName, Paralog, PDDict[Locus][Paralog][SeqName]])
			OutList.append(Line+"\n")
OutFileName = OutFolder+OutFilePre+'pardict.txt'
OutFileWriting(OutFileName, OutList)


#writing the scripts to move the files and to align the new paralog files:
#***Note that this script does not actually seem to align the new paralog files.***
ScriptList1 = ['#! /bin/bash\n\n']
for Locus in RenamedLocusList:
	if (InFolder == OutFolder) and (InFilePre == OutFilePre):
		Line = "mv "+InFolder+InFilePre+Locus+"_allbest_al.fa "+InFolder+InFilePre+Locus+"_allbest_al_old.fa\n"
		Line += "mv "+InFolder+InFilePre+Locus+"_allseqs_al.fa "+InFolder+InFilePre+Locus+"_allseqs_al_old.fa\n"
		Line += "mv "+InFolder+"RAxML_bipartitions."+InFilePre+"o150_"+Locus+" "+InFolder+"RAxML_bipartitions."+InFilePre+"o150_"+Locus+"_old\n"
		Line += "mv "+InFolder+"RAxML_bipartitions."+InFilePre+"ab_"+Locus+" "+InFolder+"RAxML_bipartitions."+InFilePre+"ab_"+Locus+"_old\n"
		ScriptList1.append(Line)
		Line = ScriptPath+"treerenamer.py "+InFolder+InFilePre+Locus+"_allbest_al_old.fa "+OutFolder+OutFilePre+Locus+"_allbest_al.fa "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+InFilePre+Locus+"_allseqs_al_old.fa "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+InFilePre+Locus+"_allover150_al_old.fa "+OutFolder+OutFilePre+Locus+"_allover150_al.fa "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+"RAxML_bipartitions."+InFilePre+"o150_"+Locus+"_old "+OutFolder+"RAxML_bipartitions."+OutFilePre+"o150_"+Locus+" "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+"RAxML_bipartitions."+InFilePre+"ab_"+Locus+"_old "+OutFolder+"RAxML_bipartitions."+OutFilePre+"ab_"+Locus+" "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		ScriptList1.append(Line)
	else:
		Line = ScriptPath+"treerenamer.py "+InFolder+InFilePre+Locus+"_allbest_al.fa "+OutFolder+OutFilePre+Locus+"_allbest_al.fa "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+InFilePre+Locus+"_allseqs_al.fa "+OutFolder+OutFilePre+Locus+"_allseqs_al.fa "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+InFilePre+Locus+"_allover150_al.fa "+OutFolder+OutFilePre+Locus+"_allover150_al.fa "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+"RAxML_bipartitions."+InFilePre+"o150_"+Locus+" "+OutFolder+"RAxML_bipartitions."+OutFilePre+"o150_"+Locus+" "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		Line += ScriptPath+"treerenamer.py "+InFolder+"RAxML_bipartitions."+InFilePre+"ab_"+Locus+" "+OutFolder+"RAxML_bipartitions."+OutFilePre+"ab_"+Locus+" "+OutFolder+OutFilePre+Locus+"_renaming_key.txt\n"
		ScriptList1.append(Line)
OutFileName1 = OutFolder+OutFilePre+"seq_renaming_script.sh"
OutFileWriting(OutFileName1, ScriptList1)
