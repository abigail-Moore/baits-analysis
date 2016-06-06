#! /usr/bin/env python

#tassembly_to_loci.py version 1.0 1 June 2015 Abby Moore
#This script parses the output from running an assembler on the sequences that
#blasted to individual loci (using a shell script written by tblast_to_masurca.py/etc.)
#to make a fasta file that contains all of the contigs for each locus (from all 
#individuals).
#The LocusFile (prefix_Locus_List.txt written by tblast_to_masurca/etc.) should have two columns:
#Locus [0]
#comma-separated list of individuals that have sequences for that locus [1]

import sys#We want to be able to talk to the command line
import re#We need regular expressions
from collections import defaultdict #I am sure that these will come in handy somewhere
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from numpy import mean#to be able to calculate averages

#Example:
'''
tassembly_to_loci.py InFolder SeqFilePre LocusFileName OutFolder OutFilePre BlastFolder BlastDBPost BlastFilePre AMode [spades, masurca, minimo, ssake]
'''

print ("%s\n" % (" ".join(sys.argv)))

Usage ='''
tassembly_to_loci.py version 1.0
This is a script to combine the assembled contigs from multiple individuals for a
single locus so they can be analyzed further.
tassembly_to_loci.py [folder with the folders of sequences arranged according to locus]
[prefix for sequence files--only for minimo, or "none" if none]
[file with the list of loci] [output folder]
[output file prefix, or none if no prefix] [folder of blast databases] [postscript
for blast database file names] [prefix for blast output] [assembly mode: one of
spades, masurca, minimo, or ssake]'''

AModeList = ['spades', 'masurca', 'minimo', 'ssake']

if len(sys.argv) != 10:
	sys.exit(Usage)
else:
	InFolder = sys.argv[1]
	SeqFilePre = sys.argv[2]
	if SeqFilePre == "none":
		SeqFilePre = ""
	LocusFileName = sys.argv[3]
	OutFolder = sys.argv[4]
	if sys.argv[5] == "none":
		OutFilePre = ""
	else:
		OutFilePre = sys.argv[5]
	BlastFolder = sys.argv[6]
	BlastDBPost = sys.argv[7]
	if BlastDBPost == "none":
		BlastDBPost = ""
	BlastFilePre = sys.argv[8]
	AMode = sys.argv[9]
	if (AMode in AModeList) == False:
		sys.exit("You wanted %s as the assembly mode, but the assembly mode can only be one of the following: %s." % (AMode, " ".join(AModeList)))


#####################################################################################
#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
#from tbaits_intron_removal.py
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	#print("Output file %s written.\n" % (FileName))
	#sys.stderr.write("Output file %s written.\n" % (FileName))

######################################################################################

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

#various dictionaries to fill out:
LTInfoDict = defaultdict(list)#LTInfoDict[Locus] = [list of contig lengths]
LIInfoDict = defaultdict(dict)#LIInfoDict[Locus][Ind] = [list of contig lengths]

#getting the sequences for that locus from the various contigs.fasta files
for Locus in LocusDict:
	LIInfoDict[Locus] = defaultdict(list)
	ContigDict[Locus] = defaultdict(dict)
	NumSeqs = 0
	for Ind in LocusDict[Locus]:
		InFileExists = "no"
		#see if the a file exists for that individual for that locus for that format:
		if AMode == "masurca":
			InFileName = InFolder+Locus+"/"+Locus+"_"+Ind+"/CA/10-gapclose/genome.ctg.fasta"
			try:
				InFile = open(InFileName, 'rU')
				InFile.close()
				InFileExists = "yes"
				ContigPre = "gc"
			except IOError:
				InFileName = InFolder+Locus+"/"+Locus+"_"+Ind+"/CA/9-terminator/genome.ctg.fasta"
				try:
					InFile = open(InFileName, 'rU')
					InFile.close()
					InFileExists = "yes"
					ContigPre = "ter"
				except IOError:
					"sadly, there are no contigs of any sort for this individual"
		elif (AMode == "spades"):
			InFileName = InFolder+Locus+"/"+Ind+"/contigs.fasta"
			try:
				InFile = open(InFileName, 'rU')
				InFile.close()
				InFileExists = "yes"
				ContigPre = "final"
			except IOError:
				InFileName = InFolder+Locus+"/"+Ind+"/K55/final_contigs.fasta"
				try:
					InFile = open(InFileName, 'rU')
					InFile.close()
					InFileExists = "yes"
					ContigPre = "K55"
				except IOError:
					InFileName = InFolder+Locus+"/"+Ind+"/K33/final_contigs.fasta"
					try:
						InFile = open(InFileName, 'rU')
						InFile.close()
						InFileExists = "yes"
						ContigPre = "K33"
					except IOError:
						InFileName = InFolder+Locus+"/"+Ind+"/K21/final_contigs.fasta"
						try:
							InFile = open(InFileName, 'rU')
							InFile.close()
							InFileExists = "yes"
							ContigPre = "K21"
						except IOError:
							"sadly, there are no contigs of any sort for this individual"
		elif AMode == "minimo":
			InFileName = InFolder+SeqFilePre+Locus+"_"+Ind+"_seqs-contigs.fa"
			try:
				InFile = open(InFileName, 'rU')
				InFile.close()
				InFileExists = "yes"
			except IOError:
				"sadly, there are no contigs of any sort for this individual"
		elif AMode == "ssake":
			InFileName = InFolder+SeqFilePre+Locus+"_"+Ind+"ssake.out.contigs"
			try:
				InFile = open(InFileName, 'rU')
				InFile.close()
				InFileExists = "yes"
			except IOError:
				"sadly, there are no contigs of any sort for this individual"
		if InFileExists == "yes":
			for record in SeqIO.parse(InFileName, 'fasta'):
				#only take sequences of a certain length, but I don't know what this should be.  The important thing is to avoid the unassembled contigs.
				SeqLen = len(str(record.seq))
				if SeqLen > 200:
					LTInfoDict[Locus].append(SeqLen)
					LIInfoDict[Locus][Ind].append(SeqLen)
					if (AMode == "spades") or (AMode == "masurca"):
						SeqName = Ind+"-"+ContigPre+"-"+record.id
					elif AMode == "minimo":
						#note that record.id is only 7!!
						#record.description: 7 len=99 nreads=2 cov=1.89899
						SeqNameTemp = record.description.split(" ")
						SeqName = Ind+"-Min-"+"_".join(SeqNameTemp)
						MinRe = re.compile("=")
						SeqName = MinRe.sub("-",SeqName) 
					elif AMode == "ssake":
						#record.id: contig1|size160|read17|cov9.20
						SeqName = Ind+"-ssake-"+record.id
					ContigDict[Locus][Ind][SeqName] = str(record.seq)
					NumSeqs += 1
	ContigNumDict[Locus] = NumSeqs
	print("%d contigs were found in a total of %d individuals for locus %s.\n" % (NumSeqs, len(ContigDict[Locus].keys()), Locus))
	sys.stderr.write("%d contigs were found in a total of %d individuals for locus %s.\n" % (NumSeqs, len(ContigDict[Locus].keys()), Locus))

#################################################################################

#writing some output files with information about the loci
#LTInfoDict[Locus] = [list of contig lengths]
#LIInfoDict[Locus][Ind] = [list of contig lengths]
LocusList = sorted(LTInfoDict.keys())
#contigs/locus (total, for the entire group), individuals/locus, mean contig length per locus
LocusInfoList = ['Locus\tNumber_of_Contigs\tNumber_of_Individuals\tMean_Contig_Length\n']
for Locus in LocusList:
	Line = Locus+"\t"+str(len(LTInfoDict[Locus]))+"\t"+str(len(LIInfoDict[Locus].keys()))+"\t"+str(mean(LTInfoDict[Locus]))+"\n"
	LocusInfoList.append(Line)
OutFileName = OutFolder+OutFilePre+"Locus_Contig_Info.txt"
OutFileWriting(OutFileName, LocusInfoList)

IndList = [ ]
for Locus in LIInfoDict:
	IndList += LIInfoDict[Locus].keys()
IndList = sorted(list(set(IndList)))
#contigs/individual/locus
NumIndContigsList = ['Locus\t'+'\t'.join(IndList)+'\n']
#mean contig length per individual per locus
IndContigLenList = ['Locus\t'+'\t'.join(IndList)+'\n']
#list of all contigs lengths per individual per locus
AllIndContigsList  = ['Locus\t'+'\t'.join(IndList)+'\n']
for Locus in LIInfoDict:
	LineNI = Locus
	LineIL = Locus
	LineAI = Locus
	for Ind in IndList:
		if LIInfoDict[Locus][Ind] == [ ]:
			LineNI += "\t"
			LineIL += "\t"
			LineAI += "\t"
		else:
			LineNI += "\t"+str(len(LIInfoDict[Locus][Ind]))
			LineIL += "\t"+str(mean(LIInfoDict[Locus][Ind]))
			LineAI += "\t"+",".join(str(CLen) for CLen in LIInfoDict[Locus][Ind])
	LineNI += '\n'
	LineIL += '\n'
	LineAI += '\n'
	NumIndContigsList.append(LineNI)
	IndContigLenList.append(LineIL)
	AllIndContigsList.append(LineAI)
OutFileName = OutFolder+OutFilePre+"Num_Contigs_per_Individual.txt"
OutFileWriting(OutFileName, NumIndContigsList)
OutFileName = OutFolder+OutFilePre+"Mean_Contig_Len_per_Ind.txt"
OutFileWriting(OutFileName, IndContigLenList)
OutFileName = OutFolder+OutFilePre+"All_Contigs_per_Ind.txt"
OutFileWriting(OutFileName, AllIndContigsList)

#writing everything to the script
BlastList1 = ['#! /bin/bash\n']
BlastList2 = [ ]
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
		BlastList2.append(Line)
		Line = BlastFilePre+Locus+".out\t"+OutFilePre+Locus+".fa\n"
		BlastFileList.append(Line)
print("Contigs for %d loci were written to files with names such as %s.\n" % (len(ContigDict.keys()), OutFileName))
sys.stderr.write("Contigs for %d loci were written to files with names such as %s.\n" % (len(ContigDict.keys()), OutFileName))

OutFileName1 = OutFolder+BlastFilePre+"BlastScript1.sh"
OutFileName2 = OutFolder+BlastFilePre+"BlastScript2.sh"
Line = "cat "+OutFileName2+" | parallel --joblog "+OutFolder+BlastFilePre+"parallel_log.log\n"
BlastList1.append(Line)

OutFileWriting(OutFileName1, BlastList1)
OutFileWriting(OutFileName2, BlastList2)
print("The script for blasting these files against the exon databases was written to %s.\n" % (OutFileName1))
sys.stderr.write("The script for blasting these files against the exon databases was written to %s.\n" % (OutFileName1))

OutFileName = OutFolder+BlastFilePre+"BlastFileList.txt"
OutFileWriting(OutFileName, BlastFileList)
print("The list of blast output files and their corresponding fasta files was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of blast output files and their corresponding fasta files was written to %s.\n" % (OutFileName))
