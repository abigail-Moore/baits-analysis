#! /usr/bin/env python

#tundivcontigs_combiner.py
#This script just takes all of the undivcontig fasta files and the Loci_to_Redo info files
#produced by tcontigs_to_fixed_paralogs.py when it is run in final mode and combines them.
#Then it writes a script to align the contigs to the backbone alignment and they are ready
#for analysis by tcontig_selection.py

#It reads the list of paralogs from the OutFilePre_Loci_to_Redo.txt file, produced by tcontigs_to_fixed_paralogs.py
'''
(no header, tab-delimitted)
amk: Locus [0]
amk1: Paralog [1]
'''
#It also needs a Groups_List file:
'''
(no header, tab-delimitted)
Anacampserotaceae: Group [0]
Anacs/none: GroupName (or "none") [1]
'''

'''
tundivcontigs_combiner.py GFileName SeqFilePath SeqFolderPre SeqFolderPost SeqFilePre AlFolder AlFilePre AlFilePost OutFolder OutFilePre NCores
tundivcontigs_combiner.py ~/transcriptomes/general/Group_List_sand2.txt ~/transcriptomes/ none spades_amk_sequences/ amks2_ ~/transcriptomes/sandbox/amk/amk_round1/ amkrnd1_ _allseqs_al ~/transcriptomes/sandbox/amk/amk_contigsplit/ amkcs_ 2
'''

import sys
from collections import defaultdict
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects


Usage = '''
tundivcontigs_combiner.py takes all of the undivided contigs from a given locus and
combines them for alignment and subsequent analysis by tcontig_selection.py
[text file with the names of the different folders, with the name of the sequences
folder in the first column and the shortened name for the subfolders in the 2nd]
[directory in which all of these files are found]
[anything between the name of the main folder and the name of the subfolder or 
"none", if nothing]
[anything after the name of the subfolder, or "none", if nothing--this can just
be "/", if there is no suffix, but there is still a folder]
[prefix for the sequence files, or "none" if none]
[folder of the backbone alignment files]
[prefix of the backbone alignment files, or "none" if none]
[suffix of the backbone alignment files, or "none" if none]
[output folder]
[prefix for output files, or "none" if none]
[number of cores for parallelization]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 12:
	sys.exit("ERROR!  This script requires 11 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
GFileName = sys.argv[1]
SeqFilePath = sys.argv[2]
if SeqFilePath[-1] != "/":
	SeqFilePath += "/"
SeqFolderPre = sys.argv[3]
if SeqFolderPre == "none":
	SeqFolderPre = ""
SeqFolderPost = sys.argv[4]
if SeqFolderPost[-1] != "/":
	SeqFolderPost += "/"
if SeqFolderPost == "none/":
	SeqFolderPost = ""
SeqFilePre = sys.argv[5]
if SeqFilePre == "none":
	SeqFilePre = ""
AlFolder = sys.argv[6]
if AlFolder[-1] != "/":
	AlFolder += "/"
AlFilePre = sys.argv[7]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[8]
if AlFilePost == "none":
	AlFilePost = ""
OutFolder = sys.argv[9]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[10]
if OutFilePre == "none":
	OutFilePre = ""
NCores = sys.argv[11]

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
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is GroupDict

#LocusGroupCombiner takes a dictionaries of the form SDict[Locus][Paralog][Group] = list of ParalogTemp and GDict[Group] = GAbbrev
#and makes separate files for all of the sequences in each paralog.  It also makes a list of all of the
#sequence names for each group and paralog.
#modified from LocusSeqGetter in tparalog_combiner.py
def LocusGroupCombiner(SDict,GDict,FilePath,FolderPre,FolderPost,FilePre,FilePost,SeqFormatIn,OFolder,OFilePre,OFilePost,SeqFormatOut):
	NumInFiles = 0
	NumOutFiles = 0
	DictTemp = defaultdict(dict)#DictTemp[Locus][Paralog][Group] = list of contigs
	for Locus in SDict:
		for Paralog in SDict[Locus]:
			OutFileName = OFolder+OFilePre+Paralog+OFilePost
			OutFile = open(OutFileName, 'w')
			NumOutFiles += 1
			DictTemp[Locus][Paralog] = defaultdict(list)
			for Group in SDict[Locus][Paralog]:
				#DictTemp[Locus][Paralog][Group] = [ ]
				for ParalogTemp in SDict[Locus][Paralog][Group]:
					InFileName = FilePath+Group+"/"+FolderPre+GDict[Group]+"/"+FolderPost+FilePre+ParalogTemp+FilePost
					if (FolderPre == "") and (GDict[Group] == ""):
						InFileName = FilePath+Group+"/"+FolderPost+FilePre+ParalogTemp+FilePost
					InFile = open(InFileName, 'rU')
					for record in SeqIO.parse(InFile, SeqFormatIn):
						SeqName = record.id+"_"+ParalogTemp
						Record1 = SeqRecord(seq=Seq(str(record.seq), IUPAC), id = SeqName, description = "")
						SeqIO.write(Record1, OutFile, SeqFormatOut)
						DictTemp[Locus][Paralog][Group].append(SeqName)
					InFile.close()
					print("File %s for locus %s was read.\n" % (InFileName, Locus))
					NumInFiles += 1
			OutFile.close()
	print("%d sequence files for %d loci were read, with names such as %s.\n" % (NumInFiles, len(SDict), InFileName))
	sys.stderr.write("%d sequence files for %d loci were read, with names such as %s.\n" % (NumInFiles, len(SDict), InFileName))
	print("%d combined sequence files were written, with names such as %s.\n" % (NumOutFiles, OutFileName))
	sys.stderr.write("%d combined sequence files were written, with names such as %s.\n" % (NumOutFiles, OutFileName))
	return DictTemp
	#DictTemp becomes ContigNameDict

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
	
	
#####################################################################################################################################

#make the dictionary of the taxonomic groups and their shortened names
#GroupDict[Group] = GroupAbb
GroupDict = DictFromFile(GFileName)
for Group in GroupDict:
	if GroupDict[Group] == "none":
		GroupDict[Group] = ""

#Now reading the Loci_to_Redo text files and combining them into RedoDict
RedoDict = defaultdict(dict)#RedoDict[Locus][Paralog][Group] = list of subparalogs
for Group in GroupDict:
	InFileName = SeqFilePath+Group+"/"+SeqFolderPre+GroupDict[Group]+"/"+SeqFolderPost+SeqFilePre+"Loci_to_Redo.txt"
	if (SeqFolderPre == "") and (GroupDict[Group] == ""):
		InFileName = SeqFilePath+Group+"/"+SeqFolderPost+SeqFilePre+"Loci_to_Redo.txt"
	try:
		InFile = open(InFileName, 'rU')
		for Line in InFile:
			Line = Line.strip('\r').strip('\n').split('\t')
			Locus = Line[0]
			ParalogTemp = Line[1]
			Paralog = ParalogTemp.split("_P")[0]
			#I am unsure how to make this dictionary, and this is a very clumsy way to do it.
			try:
				RedoDict[Locus][Paralog][Group].append(ParalogTemp)
			except KeyError:
				RedoDict[Locus][Paralog] = defaultdict(dict)
				RedoDict[Locus][Paralog][Group] = [ParalogTemp]
			except AttributeError:
				RedoDict[Locus][Paralog][Group] = [ParalogTemp]
		InFile.close()
		#print("File %s was read.\n" % (InFileName))		
	except IOError:
		"no file, do nothing"
		#print("File %s does not exist.\n" % (InFileName))

#Now reading the sequence files and combining them.

(ContigNameDict) = LocusGroupCombiner(RedoDict, GroupDict, SeqFilePath, SeqFolderPre, SeqFolderPost, "undivcontigs_"+SeqFilePre, ".fa", "fasta", OutFolder, OutFilePre, ".fa", "fasta")

#Writing the ContigNameDict (which contigs belong to which groups)
OutFileName = OutFolder+OutFilePre+"Contigs_to_Redo.txt"
OutList = [ ]
for Locus in ContigNameDict:
	for Paralog in ContigNameDict[Locus]:
		for Group in ContigNameDict[Locus][Paralog]:
			Line = Locus+"\t"+Paralog+"\t"+Group+"\t"+",".join(ContigNameDict[Locus][Paralog][Group])+"\n"
			OutList.append(Line)
OutFileWriting(OutFileName, OutList)

#Writing the script to align the data.
#the main script
OutFileName1 = OutFolder+OutFilePre+"Redo_Contigs_Alignment.sh"
#the subscript to be called by parallel
OutFileName2 = OutFolder+OutFilePre+"Redo_Alignment_Subscript.sh"
OutList1 = ["#! /bin/bash\n"]
OutList2 = ["#! /bin/bash\n"]
for Locus in RedoDict:
	for Paralog in RedoDict[Locus]:
		Line = "rm "+OutFolder+OutFilePre+Paralog+"_allseqs_al.fa\n"
		OutList1.append(Line)
		Line = "mafft --localpair --addfragments "+OutFolder+OutFilePre+Paralog+".fa --quiet --thread -1 "+AlFolder+AlFilePre+Locus+AlFilePost+".fa > "+OutFolder+OutFilePre+Paralog+"_allseqs_al.fa\n"
		OutList2.append(Line)
Line = "chmod u+x "+OutFileName2+"\n"
Line += "cat "+OutFileName2+" | parallel --jobs "+NCores+" --progress\n"
OutList1.append(Line)
OutFileWriting(OutFileName1, OutList1)
OutFileWriting(OutFileName2, OutList2)
