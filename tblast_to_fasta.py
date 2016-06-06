#! /usr/bin/env python

#tblast_to_fasta.py version 1.0 25 March 2015
#This script reads the Seqs_to_Loci.txt output 
#from tbaits_blastn_parse.py or tbaits_blast_parse.py
#and writes files for each locus for each individual.
#It expects file with the header:
#Locus_Name	Individual_Name	Number_of_Hits	Sequence_Names
#and subsequent lines like (tab-delimited):
#18S	Lewisia_cotyledon_cotyledon_14	1	Lewisia_cotyledon_cotyledon_14_1308_1249_14445_R1;Lewisia_cotyledon_cotyledon_14_1311_17462_24991_R1....
#[0]: locus name
#[1]: individual name
#[2]: number of BLAST hits these sequences had
#[3]: semi-colon-delimited list of sequence names (from fasta file)
#These sequence names will need to be converted to the actual sequence names in the BLAST file to find the sequences

import sys #to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels
import gzip #We want to be able to open zipped files.
from itertools import izip #We want to be able to look at two files simultaneously (without having to
#have them both fully in memory at the same time)

Usage ='''
tblast_to_fasta.py version 1.0 makes new fasta files for individual loci from
the blast output parsed by tbaits_blastn_parse.py or tbaits_blast_parse.py
tblast_to_fasta.py [infile name] [sequence folder] [prefix for sequence files]
[out folder] [out file prefix or "none" for none]
[mode: "together" if all reads for the same locus should be analyzed together, 
even if they come from different individuals or "separate" if reads for each locus
should be analyzed separately for the different individuals]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 7:
	sys.exit(Usage)
else:
	InFileName = sys.argv[1]
	SeqFolder = sys.argv[2]
	SeqFilePre = sys.argv[3]
	OutFolder = sys.argv[4]
	if sys.argv[5] == "none":
		OutFilePre = ""
	else:
		OutFilePre = sys.argv[5]
	OutMode = sys.argv[6]
	if (OutMode != "together") and (OutMode != "separate"):
		sys.exit("The output mode can only be 'together' or 'separate', but you wrote %s.\n" % (OutMode))

if SeqFolder[-1] != "/":
	SeqFolder += "/"

SeqNameDict = defaultdict(dict) #The dictionary that will have the sequence names from each locus
OutScriptDict = defaultdict(dict) #The dictionary that will be used to make the outfile.

#read input file and make defaultdict of form SeqNameDict[Ind][SeqName] = Locus
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	if Line[0] != "Locus_Name":
		Locus = Line[0]
		Ind = Line[1]
		SeqList = Line[3].split(';')
		for SeqName in SeqList:
			SeqNameDict[Ind][SeqName] = Locus
InFile.close()

IndList = SeqNameDict.keys()
IndList = sorted(IndList)

print("Information on sequences from %d individuals was read from the file %s.\n" % (len(IndList), InFileName))
sys.stderr.write("Information on sequences from %d individuals was read from the file %s.\n" % (len(IndList), InFileName))

#for each individual at once:
#read fastq files and make defaultdict indseqs[seqname]['R1' or 'R2'] = fastq sequence
#while doing this, check both sequences to see if they are entirely Ns, in which case
#write "error" for the fastq sequence
for Ind in IndList:
	IndSeqs = defaultdict(dict) #dictionary of sequences for that individual
	IndSeqs['R1'] = defaultdict(dict)
	IndSeqs['R2'] = defaultdict(dict)
	IndSeqs['S'] = defaultdict(dict)
	SeqName = ""
	BothGood = 0
	OneGood = 0
	BothBad = 0
	InFileName1 = SeqFolder+SeqFilePre+Ind+"_R1.fa"
	InFile1 = open(InFileName1, 'rU')
	InFileName2 = SeqFolder+SeqFilePre+Ind+"_R2.fa"
	InFile2 = open(InFileName2, 'rU')
	for Line1, Line2 in izip(InFile1, InFile2):
		Line1 = Line1.strip('\n').strip('\r')
		Line2 = Line2.strip('\n').strip('\r')
		if Line1[0] == ">":
			SeqNameNew = Line1[1:-3]
			try:
				Locus = SeqNameDict[Ind][SeqName]
				#check to see if we want the sequence before adding it to the dictionary of sequences to be written
				NumNs1 = Line1.count('N')
				NumNs2 = Line2.count('N')
				if (NumNs1 < 3) and (NumNs2 < 3):
					IndSeqs['R1'][Locus][SeqName] = Seq1
					IndSeqs['R2'][Locus][SeqName] = Seq2
					BothGood += 1
				elif (NumNs1 < 3) and (NumNs2 > 3):
					IndSeqs['S'][Locus][SeqName] = Seq1
					OneGood += 1
				elif (NumNs1 > 3) and (NumNs2 < 3):
					IndSeqs['S'][Locus][SeqName] = Seq2
					OneGood += 1
				else:
					BothBad += 1
			except KeyError:
				"do nothing"
			SeqName = SeqNameNew
			Seq1 = ""
			Seq2 = ""
		else:
			Seq1 += Line1
			Seq2 += Line2
	InFile1.close()
	InFile2.close()			
	print("%s had %d sequences where both reads were good,\n" % (Ind, BothGood))
	print("%d sequences where one read was good, and %d sequences where\n" % (OneGood, BothBad))
	print("neither read was good, out of %d total blast hits.\n" % (len(SeqNameDict[Ind].keys())))
	sys.stderr.write("%s had %d sequences where both reads were good,\n" % (Ind, BothGood))
	sys.stderr.write("%d sequences where one read was good, and %d sequences where\n" % (OneGood, BothBad))
	sys.stderr.write("neither read was good, out of %d total blast hits.\n" % (len(SeqNameDict[Ind].keys())))
	#now to write the sequences to files
	NumFiles = 0
	for Locus in IndSeqs['R1']:
		if OutMode == "separate":
			OutFileName = OutFolder+OutFilePre+Locus+"_"+Ind+".fa" 
			OutFile = open(OutFileName, 'w')
			OutScriptDict[Locus][Ind]=OutFileName
		elif OutMode == "together":
			OutFileName = OutFolder+OutFilePre+Locus+".fa"
			OutFile = open(OutFileName, 'a')
			OutScriptDict[Locus] = OutFileNam 
		for SeqName in IndSeqs['R1'][Locus]:
			OutLine = ">"+SeqName+"_R1\n"+IndSeqs['R1'][Locus][SeqName]+"\n"
			OutFile.write(OutLine)
			OutLine = ">"+SeqName+"_R2\n"+IndSeqs['R2'][Locus][SeqName]+"\n"
			OutFile.write(OutLine)
		for SeqName in IndSeqs['S'][Locus]:
			OutLine = ">"+SeqName+"_S\n"+IndSeqs['S'][Locus][SeqName]+"\n"
			OutFile.write(OutLine)
		OutFile.close()
		NumFiles += 1
	print("Sequences were written to %d files, with names such as %s.\n" % (NumFiles, OutFileName))
	sys.stderr.write("Sequences were written to %d files, with names such as %s.\n" % (NumFiles, OutFileName))

#Making a script to analyze the sequences using Minimo:
OutList = [ ]
Line = "#! /bin/bash\n\n"
OutList.append(Line)
Line = "mkdir "+OutFolder+OutFilePre+"Minimo_Assemblies/\n"
OutList.append(Line)
for Locus in OutScriptDict.keys():
	if OutMode == "separate":
		for Ind in OutScriptDict[Locus]:
			Line = "Minimo "+OutScriptDict[Locus][Ind]+" -D FASTA_EXP=1 -D MIN_LEN=40 -D MIN_IDENT=90\n"
			OutList.append(Line)
	elif OutMode == "together":
		Line = "Minimo "+OutScriptDict[Locus]+" -D FASTA_EXP=1 -D MIN_LEN=40 -D MIN_IDENT=90\n"
		OutList.append(Line)

OutFileName = OutFolder+"Minimo_script_"+OutMode+".sh"
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

if OutMode == "separate":
	print("The shell script to analyze the individuals separately using Minimo is %s.\n" % (OutFileName))
	sys.stderr.write("The shell script to analyze the individuals separately using Minimo is %s.\n" % (OutFileName))
elif OutMode == "together":
	print("The shell script to analyze the all individuals together using Minimo is %s.\n" % (OutFileName))
	sys.stderr.write("The shell script to analyze the all individuals together using Minimo is %s.\n" % (OutFileName))

LocusList = OutScriptDict.keys()
LocusList = sorted(LocusList)

OutFileName = OutFolder+OutFilePre+"Locus_List.txt"
OutFile = open(OutFileName, 'w')
for Locus in LocusList:
	OutFile.write(Locus)
	OutFile.write("\n")
OutFile.close()

print("The list of loci was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The list of loci was written to the file %s.\n" % (OutFileName))
