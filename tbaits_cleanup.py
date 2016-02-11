#! /usr/bin/env python


#tbaits_cleanup.py version 1.0 26 June 2015 Abby Moore
#This script gets rid of the files I don't need or are very easy to recreate,
#so that the disk space will not all be filled up.
#**********IMPORTANT********This is now set so that it does not deal with the
#spades contigs, because it is assuming it is run as part of the tree-testing pipeline.

import sys
from collections import defaultdict

Usage = '''
tbaits_cleanup.py version 1.0
Removes unnecessary files from the bait analysis pathway.
tbaits_cleanup.py
[tab-delimitted file with the name of the locus in the first column and the name
of one paralog in the second column--This is produced by tparalog_combiner.py]
[text file with the names of the different folders, with the name of the sequences
folder in the first column and the shortened name for the subfolders in the 2nd]
[folder in which all data files are found]
[type of alignment contigs: spades, minimo, or mazurca]
[end of the exon file folder]
[end of the paralog file folder]
[prefix for the sequence files--not including the letter for spades/minimo/masurca
or the exons/paralogs/sequences]
[output folder]
[prefix for output files, or none, if none]
'''

#tbaits_cleanup.py ~/transcriptomes/sandbox/LPList_ssb1.txt ~/transcriptomes/logs/Group_List_ADLM.txt /users/ajm3/scratch/eedwards/ spades exons_baits paralogs_baits b ~/transcriptomes/sandbox/ trial1
#tbaits_cleanup.py LLFileName GFileName DataFolder ContigType ExonFolderPost ParFolderPost SeqFilePre OutFolder OutFilePre
ContigTypeList = ['spades', 'minimo', 'masurca']

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 10:
	sys.exit("ERROR!  This script requires 9 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
LLFileName = sys.argv[1]
GFileName = sys.argv[2]
DataFolder = sys.argv[3]
if DataFolder[-1] != "/":
	DataFolder += "/"
ContigType = sys.argv[4]
if (ContigType in ContigTypeList) == False:
	sys.exit("ERROR! You requested the contig type %s, but the only available contig types are %s." % (ConitgType, ",".join(ContigTypeList)))
ExonFolderPost = sys.argv[5]
ParFolderPost = sys.argv[6]
SeqFilePre = sys.argv[7]
OutFolder = sys.argv[8]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[9]
if OutFilePre == "none":
	OutFilePre = ""

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

#ListDictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value.  Each key has multiple values, and they
#will be made into a list
#modified from DictFromFile
def ListDictFromFile(FileName):
	TempDict = defaultdict(list)
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]].append(Line[1])
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is LPDict

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
	return

if ContigType == "spades":
	DataFolder1Pre = "s2_"
	DataFolder2Pre = "spades_"
	DataFilePre = "s"
elif ContigType == "minimo":
	DataFolder1Pre = "min_"
	DataFolder2Pre = "minimo_"
	DataFilePre = "min"
elif ContigType == "masurca":
	DataFolder1Pre = "mas_"
	DataFolder2Pre = "masurca_"
	DataFilePre = "mas"

#getting the list of loci from the file
LPDict = ListDictFromFile(LLFileName)
#make the dictionary of the taxonomic groups and their shortened names
#GroupDict[Group] = GroupAbb
GroupDict = DictFromFile(GFileName)

OutList = []
Line = "#! /bin/bash\n#SBATCH -J "+OutFilePre+"\n#SBATCH -t 8:00:00\n#SBATCH -n 1\n#SBATCH --mem=16G\n"
OutList.append(Line)
for Group in GroupDict:
	#b1_GroupDict[Group] is the blast results, which take a long time to get, so keep those
	#b2_GroupDict[Group] is the blast results, which take a long time to get, so keep those
	Line = "rm -r "+DataFolder+Group+"/s1_"+GroupDict[Group]+"\n"
	#s2_GroupDict[Group]:
	Line += "rm "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/s2_*.fastq\n"
	for Locus in LPDict:
		Line += "rm -r "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+Locus+"\n"
	#OutList.append(Line)
	Line = ""
	#spades_contigs are the concatenated contigs from the spades analyses and the results from blasting those against the blast dictionaries
	#These are used as the basis for any further analyses, so keep all of this.
	#spades_exons (often with something in the middle or at the end)
	#keep: RAxML_originalLabelledTree, RAxML_classificationLikelihoodWeights, plain .fa file and _exons_al.fa files
	for FileName in ['RAxML_classification.', 'RAxML_entropy.', 'RAxML_info.', 'RAxML_labelledTree.']:
		Line += "rm "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre+ExonFolderPost+"/"+FileName+"*\n"
	OutList.append(Line)
	Line = "rm "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre+ExonFolderPost+"/RAxML_portableTree.*\n"
	Line += "rm "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre+ExonFolderPost+"/*_exons_al.phy\n"
	Line += "rm "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre+ExonFolderPost+"/*.phy.reduced\n"
	OutList.append(Line)
	Line = ""
	#spades_paralogs
	#can probably remove the .fa but not the _pars_al.fa
	for Locus in LPDict:
		for Paralog in LPDict[Locus]:
			Line += "rm "+DataFolder+Group+"/"+DataFolder1Pre+GroupDict[Group]+"/"+DataFolder2Pre+ParFolderPost+"/"+DataFilePre+"p"+SeqFilePre+"1_"+Paralog+".fa\n"
	OutList.append(Line)
	#And obviously the corresponding things for these for the tree trial script.
	#spades_sequences
	#I need all of these, but they are just sequence files.  There is a script to make trees, but that has generally not been made.

OutFileName = OutFolder+OutFilePre+"file_removal_script.sh"
OutFileWriting(OutFileName, OutList)

