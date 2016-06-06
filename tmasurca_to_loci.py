#! /usr/bin/env python

#tmasurca_to_loci.py version 1.0 26 May 2015 Abby Moore
#This script parses the output from running masurca on the sequences that
#blasted to individual loci (using a shell script written by tblast_to_masurca.py)
#to make a fasta file that contains all of the contigs for each locus (from all 
#individuals).
#The LocusFile (prefix_Locus_List.txt written by tblast_to_masurca) should have two columns:
#Locus [0]
#comma-separated list of individuals that have sequences for that locus [1]

import sys#We want to be able to talk to the command line
from collections import defaultdict #I am sure that these will come in handy somewhere

#Example:
'''
tmasurca_to_loci.py [InFolder] [LocusFileName] [OutFolder] [OutFilePre] [BlastFolder] [BlastDBPost] [BlastFilePre]
'''

print ("%s\n" % (" ".join(sys.argv)))

Usage ='''
tmasurca_to_loci.py version 1.0
This is a script to combine the masurca contigs from multiple individuals for a
single locus so they can be analyzed further.
tmasurca_to_loci.py [folder with the folders of sequences arranged according to locus]
[file with the list of loci] [output folder]
[output file prefix, or none if no prefix] [folder of blast databases] [postscript
for blast database file names] [prefix for blast output]'''

if len(sys.argv) != 8:
	sys.exit(Usage)
else:
	InFolder = sys.argv[1]
	LocusFileName = sys.argv[2]
	OutFolder = sys.argv[3]
	if sys.argv[4] == "none":
		OutFilePre = ""
	else:
		OutFilePre = sys.argv[4]
	BlastFolder = sys.argv[5]
	BlastDBPost = sys.argv[6]
	if BlastDBPost == "none":
		BlastDBPost = ""
	BlastFilePre = sys.argv[7]

LocusDict = { } #List of loci and their associated individuals, read from LocusFileName
ContigDict = defaultdict(dict)#Dictionary of contigs arranged according to locus
#and individual: ContigDict[Locus][Ind][SeqName] = Seq
ContigNumDict = { }#The number of contigs per locus
NumSeqs = 0

if InFolder[-1:] != "/":
	InFolder += "/"
if OutFolder[-1:] != "/":
	OutFolder += "/"
if BlastFolder[-1:] != "/":
	BlastFolder += "/"

#Getting the names of the loci and the individuals that have those loci from the locus file
InFile = open(LocusFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	Locus = Line[0]
	LocusDict[Locus] = Line[1].split(",")
InFile.close()

#getting the sequences for that locus from the various contigs.fasta files
for Locus in LocusDict:
	ContigDict[Locus] = defaultdict(dict)
	NumSeqs = 0
	for Ind in LocusDict[Locus]:
		print Ind
		InFileExists = "no"
		InFileName = InFolder+Locus+"/"+Locus+"_"+Ind+"/CA/10-gapclose/genome.ctg.fasta"
		try:
			InFile = open(InFileName, 'rU')
			InFileExists = "yes"
			ContigPre = "gc"
		except IOError:
			InFileName = InFolder+Locus+"/"+Locus+"_"+Ind+"/CA/9-terminator/genome.ctg.fasta"
			try:
				InFile = open(InFileName, 'rU')
				InFileExists = "yes"
				ContigPre = "ter"
			except IOError:
				"sadly, there are no contigs of any sort for this individual"
		if InFileExists == "yes":		
			for Line in InFile:
				Line = Line.strip('\r').strip('\n')
				if Line[0] == ">":
					SeqName = Ind+"-"+ContigPre+"-"+Line[1:]
					ContigDict[Locus][Ind][SeqName] = ""
					NumSeqs += 1
				else:
					ContigDict[Locus][Ind][SeqName] += Line
			InFile.close()
	ContigNumDict[Locus] = NumSeqs
	print("%d contigs were found in a total of %d individuals for locus %s.\n" % (NumSeqs, len(ContigDict[Locus].keys()), Locus))
	sys.stderr.write("%d contigs were found in a total of %d individuals for locus %s.\n" % (NumSeqs, len(ContigDict[Locus].keys()), Locus))

#writing everything to the outfile
BlastList = ['#! /bin/bash\n']
BlastFileList = [ ]
for Locus in ContigDict.keys():
	if len(ContigDict[Locus].keys()) > 0:
		OutFileName = OutFolder+OutFilePre+Locus+".fa"
		OutList = [ ]
		NumSeqs = 0
		for Ind in ContigDict[Locus].keys():
			SeqNum = 0
			for SeqName in ContigDict[Locus][Ind].keys():
				SeqNum += 1
				NumSeqs += 1
				OutLine = ">"+SeqName+"\n"+ContigDict[Locus][Ind][SeqName]+"\n"
				OutList.append(OutLine)
		OutFile = open(OutFileName, 'w')
		for Line in OutList:
			OutFile.write(Line)
		OutFile.close()
		if NumSeqs != ContigNumDict[Locus]:
			print("ERROR!!! %d contigs were found for locus %s, but only %d of them were written to the file.\n" % \
				(NumSeqs, Locus, ContigNumDict[Locus]))
			sys.stderr.write("ERROR!!! %d contigs were found for locus %s, but only %d of them were written to the file.\n" % \
				(NumSeqs, Locus, ContigNumDict[Locus]))
		Line = "blastn -db "+BlastFolder+Locus+BlastDBPost+" -query "+OutFileName+" -out "+OutFolder+BlastFilePre+Locus+".out -outfmt 6 -task blastn\n"
		BlastList.append(Line)
		Line = BlastFilePre+Locus+".out\t"+OutFilePre+Locus+".fa\n"
		BlastFileList.append(Line)

print("Contigs for %d loci were written to files with names such as %s.\n" % (len(ContigDict.keys()), OutFileName))
sys.stderr.write("Contigs for %d loci were written to files with names such as %s.\n" % (len(ContigDict.keys()), OutFileName))

OutFileName = OutFolder+BlastFilePre+"BlastScript.sh"
OutFile = open(OutFileName,'w')
for Line in BlastList:
	OutFile.write(Line)
OutFile.close()
print("The script for blasting these files against the exon databases was written to %s.\n" % (OutFileName))
sys.stderr.write("The script for blasting these files against the exon databases was written to %s.\n" % (OutFileName))

OutFileName = OutFolder+BlastFilePre+"BlastFileList.txt"
OutFile = open(OutFileName,'w')
for Line in BlastFileList:
	OutFile.write(Line)
OutFile.close()
print("The list of blast output files and their corresponding fasta files was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of blast output files and their corresponding fasta files was written to %s.\n" % (OutFileName))


