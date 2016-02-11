#! /usr/bin/env python

#tcontigs_to_paralogs_fixed.py version 2.0 4 Nov. 2015 Abby Moore
#This script reads the alignments of the contigs that belong to various paralogs
#(made by tclade_finder.py and MAFFT) and tries to make consensus sequences
#for each individual.
#Using the SeqFolder/SeqFilePre_Contig_Groups.txt file made by tseq_placer.py
#in the following format (tab delimitted, with header row):
'''Locus [0]
Contig [1]
Paralog [2]
'''
#Version 2.0 has three modes: intermediate, final, and finalpoly

from collections import defaultdict
from collections import Counter
import sys
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

#Example:
'''
tcontigs_to_fixed_paralogs.py ParInfoFileName SeqFolder SeqFilePre OutFolder OutFilePre OutCFolder OutCFilePre AlFolder AlFilePre FilePath OGFileName Mode[intermediate, final, finalpoly] PolyListName PolyAmbigs OldFilePre
tcontigs_to_fixed_paralogs.py ~/transcriptomes/Anacampserotaceae/spades_amk_paralogs/amkp1_Contig_Groups.txt ~/transcriptomes/Anacampserotaceae/spades_amk_paralogs/ amkp1_ ~/transcriptomes/Anacampserotaceae/spades_amk_sequences/ amks1_ ~/transcriptomes/Anacampserotaceae/spades_amk_exons/ amke2_ ~/transcriptomes/sandbox/amk/amk_round1/ amkrnd1_ none ~/transcriptomes/general/outgroup_list_new.txt intermediate none 0
These two still need to be changed to conform to the current version of this script:
tcontigs_to_fixed_paralogs.py ~/transcriptomes/Montiaceae/spades_amk_paralogs/amkp2_Contig_Groups.txt ~/transcriptomes/Montiaceae/spades_amk_paralogs/ amkp2_ ~/transcriptomes/Montiaceae/spades_amk_sequences/ amks2_ same none ~/transcriptomes/sandbox/amk/amk_round2/ amkrnd2_ none ~/transcriptomes/general/outgroup_list_new.txt final none 0
tcontigs_to_fixed_paralogs.py ~/transcriptomes/Didiereaceae/spades_amk_paralogs/amkp2_Contig_Groups.txt ~/transcriptomes/Didiereaceae/spades_amk_paralogs/ amkp2_ ~/transcriptomes/Didiereaceae/spades_amk_sequences/ amks2_ same none ~/transcriptomes/sandbox/amk/amk_round2/ amkrnd2_ none ~/transcriptomes/general/outgroup_list_new.txt finalpoly ~/transcriptomes/Didiereaceae/Did_poly.txt 2.5
'''

Usage = '''
tcontigs_to_fixed_paralogs.py version 2.0 tries to assemble the contigs that 
have been classified into paralogs into sequences.  If it can't do that, if it
is in intermediate mode, it writes a script to align the contigs to the new 
backbone alignment.  If it is in final mode, it outputs the contigs according
to paralog and they can be further analyzed by tcontig_selection.py.
tcontigs_to_paralogs.py
[file with information about the paralogs]
[folder in which the sequences are found]
[prefix for the sequence files, or "none", if none]
[output folder or "same" if the same as the input folder]
[prefix for the output files or "none" if none]
[folder to put the contigs that need to be reanalyzed, or "same" if the same as 
the input folder]
[prefix for the files of contigs that need to be reanalyzed, or "none", if none]
[for phylogenetic analysis: folder containing the template alignments]
[template file prefix or "none"]
[path to the fasta_to_phylip.py script or "none" if it is in the default path]
[file with the paralog/locus (tab) outgroup name]
[mode: "intermediate", if some paralogs will be tried again with a better 
backbone tree; "final", if the best contig should be taken for each locus, if
good sequences cannot be obtained; and "finalpoly", if some individuals are 
polyploid and the error-acceptance threshold should be increased]
[list of polyploid individuals or "none", if none]
[fold difference in the number of ambiguities polyploid sequences can have, 
compared to non-polyploid sequences, or 0 if polyploid sequences are not treated
differently from polyploid sequences]
[prefix for the original contig files (from the previous round)--in case there 
are no new trees for some loci from this round]
'''

ModeList = ['intermediate', 'final', 'finalpoly']

print("%s\n" % (" ".join(sys.argv)))
sys.stderr.write("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 16:
	sys.exit("ERROR!  This script requires 15 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
ParInfoFileName = sys.argv[1]
SeqFolder = sys.argv[2]
if SeqFolder[-1] != "/":
	SeqFolder += "/"
SeqFilePre = sys.argv[3]
if SeqFilePre == "none":
	SeqFilePre = ""
OutFolder = sys.argv[4]
if OutFolder == "same":
	OutFolder = SeqFolder
elif OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[5]
if OutFilePre == "none":
	OutFilePre = ""
OutCFolder = sys.argv[6]
if OutCFolder == "same":
	OutCFolder = SeqFolder
elif OutCFolder[-1] != "/":
	OutCFolder += "/"
OutCFilePre = sys.argv[7]
if OutCFilePre == "none":
	OutCFilePre = ""
AlFolder = sys.argv[8]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[9]
if AlFilePre == "none":
	AlFilePre = ""
FilePath = sys.argv[10]
if FilePath[-1] != "/":
	FilePath += "/"
if FilePath == "none/":
	FilePath = ""
OGFileName = sys.argv[11]
Mode = sys.argv[12]
if (Mode in ModeList) == False:
	sys.exit("ERROR! The mode needs to be in the following list: %s, but you wrote %s.\n%s" % (", ".join(ModeList), Mode, Usage))
PolyListName = sys.argv[13]
if (PolyListName != "none") and ((Mode == "intermediate") or (Mode == "final")):
	print("WARNING!!  You supplied a list of polyploid individuals, but this list will not be used because the program is running in %s mode, not in finalpoly mode.\n" % (Mode))
	sys.stderr.write("WARNING!!  You supplied a list of polyploid individuals, but this list will not be used because the program is running in %s mode, not in finalpoly mode.\n" % (Mode))
if (PolyListName == "none") and (Mode == "finalpoly"):
	sys.exit("ERROR! When the script is run in finalpoly mode, a list of polyploid individuals is required.\n%s" % (Usage))
if PolyListName == "none":
	PolyList = [ ]
PolyAmbigs = float(sys.argv[14])
if (PolyAmbigs != 0) and (Mode != "finalpoly"):
	print("WARNING!!  You have a non-zero number (%f) for the number of amibiguites in polyploid sequences, but polyploid sequences will not be analyzed differently.\n" % (PolyAmbigs))

if Mode == "intermediate":
	print("The mode is %s.  This means that polyploid individuals will not be treated differently, sequences that could not be reliably classified according to paralog will not be used, \
	and a consensus sequence can have as many ambiguities as it has contigs and still be considered good.\n" % (Mode))
	sys.stderr.write("The mode is %s.  This means that polyploid individuals will not be treated differently, sequences that could not be reliably classified according to paralog will not be used, \
	and a consensus sequence can have as many ambiguities as it has contigs and still be considered good.\n" % (Mode))
elif Mode == "final":
	print("The mode is %s.  This means that polyploid individuals will not be treated differently, all sequences will be used, \
	and a consensus sequence can have twice as many ambiguities as it has contigs and still be considered good.\n" % (Mode))
	sys.stderr.write("The mode is %s.  This means that polyploid individuals will not be treated differently, all sequences will be used, \
	and a consensus sequence can have twice as many ambiguities as it has contigs and still be considered good.\n" % (Mode))
elif Mode == "finalpoly":
	print("The mode is %s.  This means that polyploid individuals can have %.2f times more ambiguities than non-polyploid individuals, all sequences will be used, \
	and a consensus sequence can have twice as many ambiguities as it has contigs and still be considered good.\n" % (Mode, PolyAmbigs))
	sys.stderr.write("The mode is %s.  This means that polyploid individuals can have %.2f times more ambiguities than non-polyploid individuals, all sequences will be used, \
	and a consensus sequence can have twice as many ambiguities as it has contigs and still be considered good.\n" % (Mode, PolyAmbigs))
OldFilePre = sys.argv[15]
	
#This makes a defaultdict with two keys, with the value being a list, from a file
#It can have a header of multiple lines (or no lines).
#The first column (column 0) is the first level key.  The second level key can be
#in any other column, and the item to be made into a list can be in any other column.
#There can be more than 3 columns in the file, but the remaining columns will be ignored.
def LevelDictList(InFileName, NumHeaderLines, Key2Column, ListColumn):
	TempDict = defaultdict(dict)
	InFile = open(InFileName, 'rU')
	LineNum = 0
	#making the dictionary from the file
	for Line in InFile:
		if LineNum >= NumHeaderLines:
			Line = Line.strip('\r').strip('\n').split('\t')
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
	#This is ParContigDict

#DictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
#from tbaits_intron_removal.py
def DictFromFile(FileName):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]] = Line[1]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	#sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is OGDict

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
	#This is PolyList

#ContigtoInd classifies the contigs belonging to a single paralog to individual.
def ContigtoInd(CDict, dsymbol):
	TempDict = defaultdict(dict)
	for Locus in CDict:
		for Paralog in CDict[Locus]:
			TempDict[Locus][Paralog] = defaultdict(list)
			for Contig in CDict[Locus][Paralog]:
				IndName = Contig.split(dsymbol)[0]
				TempDict[Locus][Paralog][IndName].append(Contig)
	return TempDict
	#This is ParIndDict

#LocusSeqGetter reads a series of sequence files and makes a dictionary of the sequences
#that have been classified according to locus.  Used first in tbaits_intron_removal.py
#ContigDict = LocusSeqGetter(SeqFileList, SeqFolder, SeqFilePre, "_allseqs_al.fa", "fasta")
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
	#sys.stderr.write("%d sequence files were read, with names such as %s.\n" % (len(FileList), FileName))
	return TempDict
	#These are the various ContigDicts

#ConSeqMaker makes a consensus sequence from a group of aligned sequences
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
			if 'a' in PosNucs:
				if 'c' in PosNucs: ConSeq += 'm'
				elif 'g' in PosNucs: ConSeq += 'r'
				elif 't' in PosNucs: ConSeq += 'w'
				else: print("ERROR!!! We have strange bases in our sequence: %s\n" % (" ".join(PosNucs)))
			elif 'c' in PosNucs:
				if 'g' in PosNucs: ConSeq += 's'
				elif 't' in PosNucs: ConSeq += 'y' 
				else: print("ERROR!!! We have strange bases in our sequence: %s\n" % (" ".join(PosNucs)))
			elif 'g' in PosNucs:
				if 't' in PosNucs: ConSeq += 'k'
				else: print("ERROR!!! We have strange bases in our sequence: %s\n" % (" ".join(PosNucs)))
			AmbigNucs += 1
			AmbigNucList.append(SeqPos)
	return (ConSeq, AmbigNucs, AmbigNucList, NSeqs, Overlap)

#AmbigContigFinder finds suspicious contigs, contigs that have more than half, or more than 10 of, the
#sites in the list of ambiguous sequence positions
def AmbigContigFinder(SeqDict, SeqList, ANList):
	ListTemp = [ ]
	for ContigName in SeqList:
		ContigSeq = SeqDict[ContigName]
		NumSites = 0
		for SeqPos in ANList:
			if ContigSeq[SeqPos] != "-":
				NumSites += 1
		if (NumSites > 0.5*len(ANList)) or (NumSites > 10):
			ListTemp.append(ContigName)
	return ListTemp

#ContigComparison compares ambiguous contigs to the remaining sequences to determine
#which one does not fit.
def ContigComparison(SeqDict, SeqNameDict, ACList, ASList):
	TempDict = defaultdict(dict)
	AlignTemp = [ ]
	for IndName in SeqNameDict:
		for ContigName in SeqNameDict[IndName]:
			SeqTemp = SeqRecord(seq=(Seq(SeqDict[ContigName])), id=ContigName),
			AlignTemp += SeqTemp
	#put them in an alignment
	AlignTemp = MultipleSeqAlignment(AlignTemp)
	#then go through each contig and determine if its ambiguous bases fit with the
	#remaining (presumably mostly good) sequences, or not
	for AContig in ACList:
		ACSeq = SeqDict[AContig]
		ContraBase = 0
		ContigPos = 0
		AmbigPosList = [ ]
		#look through each ambiguous position
		for SeqPos in ASList:
			#if that position corresponds to a base in that contig:
			if ACSeq[SeqPos] != "-":
				#add it to the list of bases in the contig
				AmbigPosList.append(SeqPos)
				#look at the other bases in the alignment for that position
				PosNucs = [ ]
				for record in AlignTemp:
					if record[SeqPos] != '-':
						PosNucs.append(record[SeqPos])
				#count the number of times each nucleotide appears
				countedPN = Counter(PosNucs)
				#make a list of "accepted" nucleotides that appear about half the time
				AcceptedNucList = [ ]
				for Nuc in countedPN:
					if countedPN[Nuc] > 0.4*len(PosNucs):
						AcceptedNucList.append(Nuc)
				#count the times the ambiguous base in the contig is a different base than the rest of the sequences
				if (ACSeq[SeqPos] in AcceptedNucList) == False:
					ContraBase += 1
					#print("#%d: Sequence Position: %d, ACNuc: %s, AcceptedNucList: %s\n" % (ContraBase, SeqPos, ACSeq[SeqPos], ",".join(AcceptedNucList)))
		TempDict[AContig]['ContraBases'] = ContraBase
		TempDict[AContig]['AmbigPosList'] = AmbigPosList
	return TempDict				

#LocusParalogSeqWriter writes sequences in the desired format to files from a
#dictionary where the top-level key is the locus, the second-level key is the group number,
#the third-level key is the sequence name, and the value is the sequence.
#The top-level key can simply be "".
#from tclade_finder.py
def LocusParalogSeqWriter(SeqDict, Folder, Prefix, Suffix, SeqFormat):
	#SeqDict[Locus][Paralog][ContigName] = Sequence
	TempDict = defaultdict(dict)#TempDict[Locus][Paralog] = file name
	for Locus in SeqDict:
		NumFiles = 0
		for Paralog in SeqDict[Locus]:
			OutFileName = Folder+Prefix+Paralog+Suffix
			TempDict[Locus][Paralog] = Prefix+Paralog
			OutFile = open(OutFileName, 'w')
			ContigNameList = sorted(SeqDict[Locus][Paralog].keys())
			for ContigName in ContigNameList:
				Record1 = SeqRecord(seq=Seq(SeqDict[Locus][Paralog][ContigName], IUPAC), id = ContigName, description = "")
				SeqIO.write(Record1, OutFile, SeqFormat)
			OutFile.close()
			NumFiles += 1
		if NumFiles > 0:
			print("%d sequence files were written for the locus %s, with names such as %s.\n" % (NumFiles, Locus, OutFileName))
			#sys.stderr.write("%d sequence files were written for the locus %s, with names such as %s.\n" % (NumFiles, Locus, OutFileName))
		else:
			print("No paralogs were found for locus %s.\n" % (Locus))
			#sys.stderr.write("No paralogs were found for locus %s.\n" % (Locus))
	return TempDict
	#This is SeqFileDict

#LocusSeqWriter, modified from LocusParalogSeqWriter
def LocusSeqWriter(SeqDict, Folder, Prefix, Suffix, SeqFormat):
	#SeqDict[Locus][ContigName] = Sequence
	TempDict = { }#TempDict[Locus] = file name
	for Locus in SeqDict:
		OutFileName = "ERROR!!!  No paralogs for this locus!!"
		OutFileName = Folder+Prefix+Locus+Suffix
		TempDict[Locus] = Prefix+Locus
		OutFile = open(OutFileName, 'w')
		for ContigName in SeqDict[Locus]:
			Record1 = SeqRecord(seq=Seq(SeqDict[Locus][ContigName], IUPAC), id = ContigName, description = "")
			SeqIO.write(Record1, OutFile, SeqFormat)
		OutFile.close()
	print("%d sequence files were written that contain contigs that need to be analyzed further, with names such as %s.\n" % (len(TempDict), OutFileName))
	#sys.stderr.write("%d sequence files were written that contain contigs that need to be analyzed further, with names such as %s.\n" % (len(TempDict), OutFileName))
	return TempDict

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

#This is based on MRScriptWriter from tbaits_introns_removal.py.  It writes a script that aligns
#the contigs that couldn't be assembled to the new alignment and places them in the new tree.
def MRScriptWriter(SeqFileDict, OutGDict, Folder, Prefix, AFolder, APre, APost, Path, Prefix2):
	OutList = ["#! /bin/bash\n"]
	for Locus in SeqFileDict:
		Line = "rm "+Folder+Prefix+Locus+"_exons_al.fa\n"
		Line += "rm "+Folder+"RAxML*"+Prefix+Locus+"\n"
		#Line += "if test -f '"+AFolder+APre+Locus+APost+".fa'; then mafft --addfragments "+Folder+Prefix+Locus+".fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+Prefix+Locus+"_exons_al.fa; else cp "+Folder+Prefix2+Locus+".fa "+Folder+Prefix+Locus+".fa; fi\n"
		Line += "mafft --addfragments "+Folder+Prefix+Locus+".fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+Prefix+Locus+"_exons_al.fa\n"
		Line += Path+"fasta_to_phylip.py "+Folder+Prefix+Locus+"_exons_al.fa\n"
		#Line += "if test -f '"+AFolder+"RAxML_bipartitions."+APre+Locus+"'; then raxmlHPC -f v -s "+Folder+Prefix+Locus+"_exons_al.phy -n "+Prefix+Locus+" -t "+AFolder+"RAxML_bipartitions."+APre+Locus+" -m GTRCAT -o "+OutGDict[Locus]+" -w "+Folder+"; else cp "+Folder+"RAxML_originalLabelledTree."+Prefix2+Locus+" "+Folder+"RAxML_originalLabelledTree."+Prefix+Locus+" && cp "+Folder+"RAxML_classificationLikelihoodWeights."+Prefix2+Locus+" "+Folder+"RAxML_classificationLikelihoodWeights."+Prefix+Locus+"; fi\n"
		Line += "raxmlHPC -f v -s "+Folder+Prefix+Locus+"_exons_al.phy -n "+Prefix+Locus+" -t "+AFolder+"RAxML_bipartitions."+APre+Locus+" -m GTRCAT -o "+OutGDict[Locus]+" -w "+Folder+"\n"
		OutList.append(Line)
	OutFileName = Folder+Prefix+"Analysis_Script.sh"
	OutFileWriting(OutFileName, OutList)
	print("The shell script for analyzing the contigs that could not be assembled was written to %s.\n" % (OutFileName))
	sys.stderr.write("The shell script for analyzing the contigs that could not be assembled was written to %s.\n" % (OutFileName))

#########Reading everything in and organizing it:

#Read the file with information about the paralogs
#ParContigDict[Locus][Paralog] = list of Contigs
ParContigDict = LevelDictList(ParInfoFileName, 1, 2, 1)
#get information about the outgroups:
OGDict = DictFromFile(OGFileName)

#get the list of polyploid individuals, if applicable
if Mode == "finalpoly":
	PolyList = CaptureColumn(PolyListName, 0)

#Make another dictionary with the contigs separated by individual
#ParIndDict[Locus][Paralog][IndName] = list of Contigs
ParIndDict = ContigtoInd(ParContigDict, "-")

#################Attempting to make good consensus sequences for each individual for each paralog:
#go through the various loci
BadContigDict = defaultdict(dict) #BadContigDict[Locus][ContigName] = ContigSeq
NumAmbigDict = defaultdict(dict)#NumAmbigDict[Locus][Paralog]['Good'/'Bad'/'Length'] = list of NumAmbigs for different ConSeqs
IndInfoDict = defaultdict(dict)#IndInfoDict[Locus][Paralog][IndName]['NumSeqs'/'NumAmbig'/'SeqLen'] = info
ContigDict = defaultdict(dict)#ContigDict[Locus][Paralog][Contig] = ContigSeq
ContigFateDict = defaultdict(dict)#ContigFateDict[Locus][Contig] = ContigFate
for Locus in ParContigDict:
	#read the sequences
	SeqFileList = [SeqFilePre+Paralog+"_pars_al.fa" for Paralog in ParContigDict[Locus].keys()]
	ContigDict[Locus] = LocusSeqGetter(SeqFileList, SeqFolder, SeqFilePre, "_pars_al.fa", "fasta")
	#ContigDict[Locus][Paralog][Contig] = ContigSeq
	for Paralog in ParIndDict[Locus]:
		IndInfoDict[Locus][Paralog] = defaultdict(dict)
		#if this is not the final stage, we do not want to analyze contigs that were not well-classified,
		#because we want to see if their classification can be improved with an improved tree.
		if (Paralog[0:5] == "Ambig") and (Mode == "intermediate"):
			#move all of the contigs for this paralog into the BadContigDict
			for IndName in ParIndDict[Locus][Paralog]:
				for Contig in ParIndDict[Locus][Paralog][IndName]:
					BadContigDict[Locus][Contig] = ContigDict[Locus][Paralog][Contig]
					ContigFateDict[Locus][Contig] = "not_using_"+Paralog
				IndInfoDict[Locus][Paralog][IndName]['NumSeqs'] = len(ParIndDict[Locus][Paralog][IndName])
				IndInfoDict[Locus][Paralog][IndName]['NumAmbig'] = 'n/a'
				IndInfoDict[Locus][Paralog][IndName]['SeqLen'] = 'n/a'
				IndInfoDict[Locus][Paralog][IndName]['Overlap'] = 'n/a'
				IndInfoDict[Locus][Paralog][IndName]['Using'] = 'no'
		else: #analyze the sequences
			NumAmbigDict[Locus][Paralog] = defaultdict(list) #the remaining bits of NumAmbigDict[Locus][Paralog] are lists
			NumAmbigDict[Locus][Paralog]['Length'] = 0 #but 'Length' is just a value
			#calculate the consensus sequence for each individual and determine how many ambiguities it has
			for IndName in ParIndDict[Locus][Paralog]:
				#print ("%s, %s, %s" % (Locus, Paralog, IndName))
				(ConSeq, NumAmbig, AmbigSitesList, NumSeqs, Overlap) = ConSeqMaker(ContigDict[Locus][Paralog], ParIndDict[Locus][Paralog][IndName])
				SeqName = IndName+"."+Paralog+"."+str(NumSeqs)+"seqs."+str(NumAmbig)+"ambig.len"+str(len(ConSeq))
				#update the length of the consensus sequence
				if NumAmbigDict[Locus][Paralog]['Length'] < len(ConSeq):
					NumAmbigDict[Locus][Paralog]['Length'] = len(ConSeq)
				#Setting the threshold for sequence acceptance:
				if Mode == "intermediate":
					SeqThreshold = NumSeqs
					#print("The most ambiguities this sequence can have is %d.\n" % (SeqThreshold))
				elif (Mode == "finalpoly") and (IndName in PolyList):
					SeqThreshold = NumSeqs*2*PolyAmbigs
					#print("Individual %s is polyploid, so it can have up to %d amibiguities.  It has %d ambiguities.\n" % (IndName, SeqThreshold, NumAmbig))
				else:
					SeqThreshold = NumSeqs*2
					#print("The most ambiguities this sequence can have is %d.\n" % (SeqThreshold))
				#There are now three possibilities of what to do with the contig:
				#1) It can be accepted as is, because it does not have many ambiguities
				if NumAmbig <= SeqThreshold:
					NumAmbigDict[Locus][Paralog]['Good'].append(NumAmbig)
					#print("Accepted!")
				#2) If it does have too many ambiguities, but it is a polyploid, then accepting it for that reason
				elif (Mode == "finalpoly") and (IndName in PolyList):
					if NumAmbig <= SeqThreshold:
						NumAmbigDict[Locus][Paralog]['Good'].append(NumAmbig)
						#print("Accepted!")
					#but considering the sequence to be bad if it has more than that
					else:
						NumAmbigDict[Locus][Paralog]['Bad'].append(NumAmbig)
						#print("Rejected!")
				#3) Or, if neither of these things are true, then looking for chimeric bits of contigs, and accepting the sequence only if it appears as
				#if the many ambiguities are due to that.
				else:
					#print("WARNING!! The sequence for individual %s at paralog %s had %d ambiguities.  It should be analyzed differently.\n" % (IndName, Paralog, NumAmbig))
					AmbigContigList = AmbigContigFinder(ContigDict[Locus][Paralog], ParIndDict[Locus][Paralog][IndName], AmbigSitesList)
					ACInfoDict = ContigComparison(ContigDict[Locus][Paralog], ParIndDict[Locus][Paralog], AmbigContigList, AmbigSitesList)
					if len(ACInfoDict.keys()) == 0:
						#print("Unfortunately, the ambiguous bases are spread throughout the sequence, so they do not appear to be due to chimeric contigs.\n")
						NumAmbigDict[Locus][Paralog]['Bad'].append(NumAmbig)
						#print("Rejected!")
					else:
						ProblemContigs = [ ]
						#First, I need to determine if there are too many problematic contigs, making the problems more likely due to multiple sequences at that locus
						#than to a couple of chimeric contigs.
						for Contig in ACInfoDict:
							if len(ACInfoDict[Contig]['AmbigPosList']) > 4:
								ProblemContigs.append(Contig)
						if len(ProblemContigs) > 0.5*NumSeqs:
							#print("Individual %s at paralog %s has %d ambiguous contigs, making it likely that it has multiple sequences for this locus.\n" % (IndName, Paralog, len(ProblemContigs)))
							NumAmbigDict[Locus][Paralog]['Bad'].append(NumAmbig)
							#print("Rejected!")
						#if it seems like the problems could be due to chimeric contigs, look through the contigs and prune out the bad ones.
						else:
							PrunedContigs = [ ]
							for Contig in ACInfoDict:
								for Contig2 in AmbigContigList:
									if Contig != Contig2:
										#determine whether the ambiguous positions in the two contigs overlap
										Contig1set = set(ACInfoDict[Contig]['AmbigPosList'])
										Contig2set = set(ACInfoDict[Contig2]['AmbigPosList'])
										IntersectSet = Contig1set.intersection(Contig2set)
										#and, if they do, if most of them overlap or if only a few overlap
										if len(IntersectSet) >= (len(Contig1set)+len(Contig2set))/4:
											#print("%s and %s share %d of the ambiguous bases.\n" % (Contig, Contig2, len(IntersectSet)))
											#when a clearly suspicuous contig is found:
											if ACInfoDict[Contig]['ContraBases'] > 3*ACInfoDict[Contig2]['ContraBases']:
												#print("The problem appears to lie with contig %s, which has %d contradictory bases, compared to %d in contig %s." % (Contig, ACInfoDict[Contig]['ContraBases'], ACInfoDict[Contig2]['ContraBases'], Contig2))
												if (Contig in PrunedContigs) == False:
													PrunedContigs.append(Contig)
													#print("The ambiguous bases from the overlapping region will be removed from that contig.\n")
													ContigSeqList = list(ContigDict[Locus][Paralog][Contig])
													NewContigTemp = [ ]
													#find the ambiguous area
													StartBase = sorted(IntersectSet)[0]
													EndBase = sorted(IntersectSet)[-1]
													for Base in range(StartBase,EndBase+1):
														#make a new sequence with just that area
														NewContigTemp.append(ContigSeqList[Base])
														#and remove the area from the current contig
														ContigSeqList[Base] = "-"
													#Then I need some type of filter to avoid empty sequences.
													#In this case, we would only accept pruned contigs that have all 4 bases, plus the gap character that we know we added.
													if len(list(set(ContigSeqList))) > 4:
														ContigDict[Locus][Paralog][Contig] = "".join(ContigSeqList)
														#For the new contigs, they need to have either all four bases or three of the four plus gap.
														if len(list(set(NewContigTemp))) > 3:
															BadContigDict[Locus][Contig+str(StartBase)+"_to_"+str(EndBase)] = "".join(NewContigTemp)
													#Any sequence we can align will fit these criteria.  Of course, it allows really short sequences, but, since we're aligning fragments at this stage, that's probably alright.
													#And if the fragments don't have a good placement in the next step, they will just be rejected, anyway.
												else:
													"do nothing"
													#print("Contig %s has already been pruned, so it will be ignored in this round.\n" % (Contig))
							#then remake the consensus sequence once the bad portions of the contigs have been removed
							(ConSeq, NumAmbig, AmbigSitesList, NumSeqs, Overlap) = ConSeqMaker(ContigDict[Locus][Paralog], ParIndDict[Locus][Paralog][IndName])
							SeqName = IndName+"."+Paralog+"."+str(NumSeqs)+"seqs."+str(NumAmbig)+"ambig.len"+str(len(ConSeq))
							#print("Individual %s now has %d ambiguities at paralog %s.\n" % (IndName, NumAmbig, Paralog))
							#try again to see if the contig passes the filter
							if NumAmbig > SeqThreshold:
								NumAmbigDict[Locus][Paralog]['Bad'].append(NumAmbig)
								#print("Rejected!")
							else:
								NumAmbigDict[Locus][Paralog]['Good'].append(NumAmbig)
								#print("Accepted!")
				IndInfoDict[Locus][Paralog][IndName]['NumSeqs'] = NumSeqs
				IndInfoDict[Locus][Paralog][IndName]['NumAmbig'] = NumAmbig
				IndInfoDict[Locus][Paralog][IndName]['SeqLen'] = len(ConSeq)
				IndInfoDict[Locus][Paralog][IndName]['Overlap'] = Overlap
				IndInfoDict[Locus][Paralog][IndName]['SeqName'] = SeqName
				IndInfoDict[Locus][Paralog][IndName]['Seq'] = ConSeq

##########Deciding on the ultimate fate of these consensus sequences:
#looking at the number of ambiguities for the various paralogs and determining which are worth analyzing further
GoodSeqDict = defaultdict(dict)#good sequences #GoodSeqDict[Locus][Paralog][SeqName] = Conseq
BadSeqDict = defaultdict(dict)#bad sequences #BadSeqDict[Locus][Paralog][SeqName] = Conseq
UndivContigDict = defaultdict(dict)#This is for the final and finalpoly modes to house the contigs that are part of rejected consensus sequences.  The contigs here are divided into paralogs.
#UndivContigDict[Locus][Paralog][Contig] = ContigSeq
for Locus in NumAmbigDict:
	GoodSeqDict[Locus] = defaultdict(dict)
	#checking each paralog in each locus
	for Paralog in NumAmbigDict[Locus]:
		#for intermediate mode, paralogs as a whole are accepted or rejected:
		if Mode == "intermediate":
			AcceptParalog = True
			#one way a paralog can be rejected is if the number of "bad" sequences is at least half as many as the number of "good" sequences
			if 2*len(NumAmbigDict[Locus][Paralog]['Bad']) > len(NumAmbigDict[Locus][Paralog]['Good']):
				AcceptParalog = False
			#the other is for the number of ambiguities to be greater than 5% of the length of the combined sequence
			elif len(NumAmbigDict[Locus][Paralog]['Bad']) > 0:
				AmbigList = NumAmbigDict[Locus][Paralog]['Bad']
				#if there are sequences with ambiguities, checking to see if there are too many, and rejecting if so
				if (20*(sorted(AmbigList)[-1]) > NumAmbigDict[Locus][Paralog]['Length']):
					AcceptParalog = False
			#if we can accept this paralog,
			if AcceptParalog == True:
				#print("%s: %s accepted!\n" % (Locus, Paralog))
				for IndName in IndInfoDict[Locus][Paralog]:
					#add the sequences to GoodSeqDict
					SeqName = IndInfoDict[Locus][Paralog][IndName]['SeqName']
					GoodSeqDict[Locus][Paralog][SeqName] = IndInfoDict[Locus][Paralog][IndName]['Seq']
					#and record that we are using the consensus sequence for each individual
					IndInfoDict[Locus][Paralog][IndName]['Using'] = 'yes'
					#recording information about the contigs:
					for Contig in ParIndDict[Locus][Paralog][IndName]:
						ContigFateDict[Locus][Contig] = "good_"+Paralog
			elif AcceptParalog == False:
				#print("%s: %s rejected!\n" % (Locus, Paralog))
				#add the contigs for that paralog back to the pile of contigs we need to analyze again
				for IndName in ParIndDict[Locus][Paralog]:
					for Contig in ParIndDict[Locus][Paralog][IndName]:
						BadContigDict[Locus][Contig] = ContigDict[Locus][Paralog][Contig]
						ContigFateDict[Locus][Contig] = "not_using_"+Paralog
					#and record that we are not using the consensus sequence for each individual
					IndInfoDict[Locus][Paralog][IndName]['Using'] = 'no'
		elif (Mode == "final") or (Mode == "finalpoly"):
			AcceptInds = [ ]
			RejectInds = [ ]
			#if all of the sequences are good, all individuals will be accepted:
			if len(NumAmbigDict[Locus][Paralog]['Bad']) == 0:
				#print("For %s: %s, all %d of the sequences are good, so they will all be accepted.\n" % (Locus, Paralog, len(NumAmbigDict[Locus][Paralog]['Good'])))
				for IndName in IndInfoDict[Locus][Paralog]:
					AcceptInds.append(IndName)
			#if more than half of the sequences are bad, all individuals will be rejected:
			elif len(NumAmbigDict[Locus][Paralog]['Bad']) > len(NumAmbigDict[Locus][Paralog]['Good']):
				print("For %s: %s, more than half (%d of %d) sequences are bad, so they will all be rejected.\n" % (Locus, Paralog, len(NumAmbigDict[Locus][Paralog]['Bad']), len(NumAmbigDict[Locus][Paralog]['Good'])+len(NumAmbigDict[Locus][Paralog]['Bad'])))
				for IndName in IndInfoDict[Locus][Paralog]:
					RejectInds.append(IndName)
			#Analyzing the sequences individually, if the situation is more complicated:
			else:
				#print("For %s: %s, the situation is more complicated." % (Locus, Paralog))
				for IndName in IndInfoDict[Locus][Paralog]:
					NumAmbig = IndInfoDict[Locus][Paralog][IndName]['NumAmbig']
					NumSeqs = IndInfoDict[Locus][Paralog][IndName]['NumSeqs']
					#deciding which individuals are rejected and accepted
					if NumAmbig <= 2*NumSeqs:
						AcceptInds.append(IndName)
						#print("%s has a good sequence, with only %d ambiguities.\n" % (IndName, NumAmbig))
					elif (IndName in PolyList) and (NumAmbig <= PolyAmbigs*2*NumSeqs):
						AcceptInds.append(IndName)
						#print("%s has a good sequence, with only %d ambiguities.\n" % (IndName, NumAmbig))
					else:
						RejectInds.append(IndName)
						#print("%s has a bad sequence, with %d ambiguities.\n" % (IndName, NumAmbig))
				print("For %s: %s, %d individuals had good sequences, and %d had bad sequences.  Perhaps these individuals are polyploid: %s?\n" % (Locus, Paralog, len(AcceptInds), len(RejectInds), ", ".join(RejectInds)))
				#recording this information in the IndInfoDict, and putting the good and bad consensus sequences in GoodSeqDict and BadSeqDict, respectively
			#print("Accepted Individuals are: %s.\n" % (" ".join(AcceptInds)))
			for IndName in AcceptInds:
				IndInfoDict[Locus][Paralog][IndName]['Using'] = 'yes'
				SeqName = IndInfoDict[Locus][Paralog][IndName]['SeqName']
				GoodSeqDict[Locus][Paralog][SeqName] = IndInfoDict[Locus][Paralog][IndName]['Seq']
				for Contig in ParIndDict[Locus][Paralog][IndName]:
					ContigFateDict[Locus][Contig] = "good_"+Paralog
			if len(RejectInds) != 0:
				#print("Rejected Individuals are: %s.\n" % (" ".join(RejectInds)))
				#add the paralog to the list of bad paralogs for this locus*******not sure we need to do this anymore
				UndivContigDict[Locus][Paralog] = defaultdict(dict)
				BadSeqDict[Locus][Paralog] = defaultdict(dict)
				for IndName in RejectInds:
					IndInfoDict[Locus][Paralog][IndName]['Using'] = 'no'
					SeqName = IndInfoDict[Locus][Paralog][IndName]['SeqName']
					BadSeqDict[Locus][Paralog][SeqName] = IndInfoDict[Locus][Paralog][IndName]['Seq']
					#and, for the bad consensus sequences, putting the contigs in a separate file for further analysis
					#add the contigs for that paralog back to the pile of contigs we need to analyze again
					for Contig in ParIndDict[Locus][Paralog][IndName]:
						UndivContigDict[Locus][Paralog][Contig] = ContigDict[Locus][Paralog][Contig]
						ContigFateDict[Locus][Contig] = "undiv_"+Paralog


##########Now we need to output all of the sequences to the appropriate files for further analyses.
#First, we want the good sequences.
SeqFileDict = LocusParalogSeqWriter(GoodSeqDict, OutFolder, OutFilePre, '.fa', 'fasta')

#Then, what we do with the bad sequences depends on the mode.
#If the mode is intermediate, then we want all of the contigs that didn't have an unambiguous home to be in one file per locus for reanalysis in the second round
#once we have a better backbone tree:
if Mode == "intermediate":
	#writing the bad contigs to a file, assuming we have bad contigs
	if len(BadContigDict.keys()) > 0:
		#writing the contigs to separate files for each locus
		BadContigFileDict = LocusSeqWriter(BadContigDict, OutCFolder, OutCFilePre, '.fa', 'fasta')
		#and making a list of those loci that need to be reanalyzed
		RedoList = [ ]
		for Locus in BadContigDict:
			Line = Locus+"\n"
			RedoList.append(Line)
		#information about the loci that need to be redone
		OutFileName = OutFolder+OutFilePre+"Loci_to_Redo.txt"
		OutFileWriting(OutFileName, RedoList)
#But, if this is the last round, we don't want the contigs for the different paralogs to be combined anymore, because
#we want to be able to go through the individual paralogs separately
if (Mode == "final") or (Mode == "finalpoly"):
	#dealing with the bad paralogs, assuming we have any
	if len(BadSeqDict.keys()) > 0:
		#So, first we want to have the "bad" consensus sequences, in case we want to look at some of them, after all
		BadSeqFileDict = LocusParalogSeqWriter(BadSeqDict, OutFolder, "badseqs_"+OutFilePre, '.fa', 'fasta')
		#And now we want separate files for each paralog for the separate contigs, so, if we can't use the sequences, we will be able to use part of them.
		UndivSeqFileDict = LocusParalogSeqWriter(UndivContigDict, OutFolder, "undivcontigs_"+OutFilePre, '.fa', 'fasta')
		#And we need a list of those files, so the next script can analyze them automatically:
		RedoList = [ ]
		for Locus in UndivSeqFileDict:
			for Paralog in UndivSeqFileDict[Locus]:
				Line = Locus+"\t"+Paralog+"\n"
				RedoList.append(Line)
		#information about the loci that need to be redone
		OutFileName = OutFolder+OutFilePre+"Loci_to_Redo.txt"
		OutFileWriting(OutFileName, RedoList)
	
##############Now to write the information files:
#a script to analyze the bad contigs (the good sequences will be combined with those from other groups by
#tparalog_combiner.py)
#but only if this is the intermediate mode
if (Mode == "intermediate") and (len(BadContigDict) > 0):
	MRScriptWriter(BadContigDict, OGDict, OutCFolder, OutCFilePre, AlFolder, AlFilePre, "_allseqs_al", FilePath, OldFilePre)

#information about the number of ambiguities per locus:
OutList = ["Locus\tParalog\t#_Ambiguities_in_Good_Seqs\t#_Ambiguities_in_Bad_Seqs\t#_Good_Seqs\t#_Bad_Seqs\n"]
OutListLP = [ ]
OutFileName = OutFolder+OutFilePre+"Paralog_Info.txt"
#OutFileNameLP = OutFolder+OutFilePre+"Locus_Paralog_List.txt"
for Locus in NumAmbigDict:
	for Paralog in NumAmbigDict[Locus]:
		Line = Locus+"\t"+Paralog+"\t"
		if len(NumAmbigDict[Locus][Paralog]['Good']) != 0: Line += ",".join(str(Num) for Num in NumAmbigDict[Locus][Paralog]['Good'])
		Line += "\t"
		if len(NumAmbigDict[Locus][Paralog]['Bad']) != 0: Line += ",".join(str(Num) for Num in NumAmbigDict[Locus][Paralog]['Bad'])
		Line += "\t"+str(len(NumAmbigDict[Locus][Paralog]['Good']))+"\t"+str(len(NumAmbigDict[Locus][Paralog]['Bad']))+"\n"
		OutList.append(Line)
		#Line = Locus+"\t"+Paralog+"\n"
		#OutListLP.append(Line)
OutFileWriting(OutFileName, OutList)
#OutFileWriting(OutFileNameLP, OutListLP)

#information about the different final sequences
OutList = ['Locus\tParalog\tIndName\t#_Sequences\tLength_of_Overlap\t#_Ambiguities\t#Seq_Length\tUsing_Seq\n']
OutFileName = OutFolder+OutFilePre+"Ind_Seq_Info.txt"
for Locus in IndInfoDict:
	for Paralog in IndInfoDict[Locus]:
		for IndName in IndInfoDict[Locus][Paralog]:
			Line = Locus+"\t"+Paralog+"\t"+IndName+"\t"+str(IndInfoDict[Locus][Paralog][IndName]['NumSeqs'])+"\t"+str(IndInfoDict[Locus][Paralog][IndName]['Overlap'])+"\t"+str(IndInfoDict[Locus][Paralog][IndName]['NumAmbig'])+"\t"+str(IndInfoDict[Locus][Paralog][IndName]['SeqLen'])+"\t"+IndInfoDict[Locus][Paralog][IndName]['Using']+"\n"
			OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#information about the contig fates
OutList = ["Locus\tContig_Name\tFate\n"]
OutFileName = OutFolder+OutFilePre+"Contig_Fates.txt"
for Locus in ContigFateDict:
	for Contig in ContigFateDict[Locus]:
		Line = Locus+"\t"+Contig+"\t"+ContigFateDict[Locus][Contig]+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)
