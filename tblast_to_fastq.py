#! /usr/bin/env python

#tblast_to_fastq.py version 1.0 10 March 2015
#This script reads the Seqs_to_Loci.txt output 
#from tbaits_blastn_parse.py or tbaits_blast_parse.py
#and writes files for each locus for each individual.
#version 1.1 23 Feb. 2016, outscript parallelized

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
#have them both fully in memory at the same time

#examples
'''
tblast_to_fastq.py InFileName SeqFolder SeqFilePre OutFolder OutFilePre OutMode[together separate]
'''

Usage ='''
tblast_to_fastq.py version 1.0 makes new fastq files for individual loci from
the blast output parsed by tbaits_blastn_parse.py or tbaits_blast_parse.py
tblast_to_fastq.py [infile name] [sequence folder] [prefix for sequence files]
[out folder] [out file prefix or "none" for none]
[mode: "together" if all reads for the same locus should be analyzed together, 
even if they come from different individuals or "separate" if reads for each locus
should be analyzed separately for the different individuals]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 7:
	sys.exit("Error!!  tblast_to_fastq.py requires 7-8 additional arguments and you supplied %d.  %s" % (len(sys.argv)-1, Usage))
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
			SeqName = SeqName.split("_")
			SeqName = "_".join(SeqName[-3:])
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
	BothGood = 0
	OneGood = 0
	BothBad = 0
	InFileName1 = SeqFolder+SeqFilePre+Ind+"_R1.fastq"
	try:
		InFile1 = open(InFileName1, 'rU')
		InFileName2 = SeqFolder+SeqFilePre+Ind+"_R2.fastq"
		InFile2 = open(InFileName2, 'rU')
	except IOError:
		InFileName1 += ".gz"
		InFile1 = gzip.open(InFileName1, 'rU')
		InFileName2 = SeqFolder+SeqFilePre+Ind+"_R2.fastq.gz"
		InFile2 = gzip.open(InFileName2, 'rU')
	LineNum = 0
	for Line1, Line2 in izip(InFile1, InFile2):
		Line1 = Line1.strip('\n').strip('\r')
		Line2 = Line2.strip('\n').strip('\r')
		SeqLine = (LineNum + 4) % 4
		#check the name to see if this is a sequence we want
		if SeqLine == 0:
			LineTemp = Line1.split(" ")[0].split(":")
			SeqName = "_".join(LineTemp[4:])
			try:
				Locus = SeqNameDict[Ind][SeqName]
				SeqWanted = True
				Seq1 = Line1 + "\n"
				Seq2 = Line2 + "\n"
			except KeyError:
				SeqWanted = False
		#check to see if the sequences are mainly Ns
		if (SeqLine == 1) and (SeqWanted == True):
			NumNs1 = Line1.count('N')
			if NumNs1 > 3:
				Seq1 = "error"
				Seq1Good = False
			else:
				Seq1 += Line1 + "\n"
				Seq1Good = True
			NumNs2 = Line2.count('N')
			if NumNs2 > 3:
				Seq2 = "error"
				Seq2Good = False
			else:
				Seq2 += Line2 + "\n"
				Seq2Good = True
		#add the third line
		elif (SeqLine == 2) and (SeqWanted == True):
			Seq1 += Line1 + "\n"
			Seq2 += Line2 + "\n"
		#add the fourth line
		elif (SeqLine == 3) and (SeqWanted == True):
			Seq1 += Line1 + "\n"
			Seq2 += Line2 + "\n"
			#check to see if we want the sequence before adding it to the dictionary of sequences to be written
			if (Seq1Good == True) and (Seq2Good == True):
				IndSeqs['R1'][Locus][SeqName] = Seq1
				IndSeqs['R2'][Locus][SeqName] = Seq2
				BothGood += 1
			elif (Seq1Good == True) and (Seq2Good == False):
				IndSeqs['S'][Locus][SeqName] = Seq1
				OneGood += 1
			elif (Seq1Good == False) and (Seq2Good == True):
				IndSeqs['S'][Locus][SeqName] = Seq2
				OneGood += 1
			else:
				BothBad += 1
		LineNum += 1
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
			OutFileName1 = OutFolder+OutFilePre+Locus+"_"+Ind+"_R1.fastq" 
			OutFile1 = open(OutFileName1, 'w')		
			OutFileName2 = OutFolder+OutFilePre+Locus+"_"+Ind+"_R2.fastq" 
			OutFile2 = open(OutFileName2, 'w')
			OutScriptDict[Locus][Ind]=[OutFileName1, OutFileName2]
		elif OutMode == "together":
			OutFileName1 = OutFolder+OutFilePre+Locus+"_R1.fastq" 
			OutFile1 = open(OutFileName1, 'a')		
			OutFileName2 = OutFolder+OutFilePre+Locus+"_R2.fastq" 
			OutFile2 = open(OutFileName2, 'a')
			OutScriptDict[Locus]["R1_R2"] = [OutFileName1, OutFileName2] 
		for SeqName in IndSeqs['R1'][Locus]:
			OutFile1.write(IndSeqs['R1'][Locus][SeqName])
			OutFile2.write(IndSeqs['R2'][Locus][SeqName])
		OutFile1.close()
		OutFile2.close()
		NumFiles += 2
	for Locus in IndSeqs['S']:
		if OutMode == "separate":
			OutFileNameS = OutFolder+OutFilePre+Locus+"_"+Ind+"_S.fastq" 
			OutFileS = open(OutFileNameS, 'w')
			try:
				OutScriptDict[Locus][Ind].append(OutFileNameS)
			except KeyError:
				OutScriptDict[Locus][Ind] = [OutFileNameS]
		elif OutMode == "together":
			OutFileNameS = OutFolder+OutFilePre+Locus+"_S.fastq" 
			OutFileS = open(OutFileNameS, 'a')
			OutScriptDict[Locus]['S'] = [OutFileNameS]
		for SeqName in IndSeqs['S'][Locus]:
			OutFileS.write(IndSeqs['S'][Locus][SeqName])
		OutFileS.close()
		NumFiles += 1
	print("Sequences were written to %d files, with names such as %s.\n" % (NumFiles, OutFileName1))
	sys.stderr.write("Sequences were written to %d files, with names such as %s.\n" % (NumFiles, OutFileName1))

#Making a script to analyze the sequences using Spades:
OutScript1 = [ "#! /bin/bash\n\n" ]
OutScript2 = [ ]
Line = "mkdir "+OutFolder+OutFilePre+"Spades_Contigs/\n"
OutScript1.append(Line)
for Locus in OutScriptDict.keys():
	Line = "mkdir "+OutFolder+Locus+"\n"
	OutScript1.append(Line)
	if OutMode == "separate":
		for Ind in IndListOut:
			OutFolderName = OutFolder+Locus+"/"+Ind
			Line = "mkdir "+OutFolderName+"\n"
			OutScript1.append(Line)
			if len(OutScriptDict[Locus][Ind]) == 2:
				Line = "spades.py -o "+OutFolderName+" --pe1-1 "+OutScriptDict[Locus][Ind][0]+" --pe1-2 "+OutScriptDict[Locus][Ind][1]+" -t 4 -m 64\n"
				OutScript2.append(Line)
			elif len(OutScriptDict[Locus][Ind]) == 3:
				Line = "spades.py -o "+OutFolderName+" --pe1-1 "+OutScriptDict[Locus][Ind][0]+" --pe1-2 "+OutScriptDict[Locus][Ind][1]+" --pe1-s "+OutScriptDict[Locus][Ind][2]+" -t 4 -m 64\n"
				OutScript2.append(Line)
			elif len(OutScriptDict[Locus][Ind]) == 1:
				Line = "spades.py -o "+OutFolderName+" --pe1-s "+OutScriptDict[Locus][Ind][0]+" -t 4 -m 64\n"
				OutScript2.append(Line)
	elif OutMode == "together":
		OutFolderName = OutFolder+Locus+"/together"
		Line = "mkdir "+OutFolderName+"\n"
		OutScript1.append(Line)
		try:
			Line = "spades.py -o "+OutFolderName+" --pe1-1 "+OutScriptDict[Locus]['R1_R2'][0]+" --pe1-2 "+OutScriptDict[Locus]['R1_R2'][1]+" --pe1-s "+OutScriptDict[Locus]['S'][0]+" -t 4 -m 64\n"
			OutScript2.append(Line)
		except KeyError:
			Line = "spades.py -o "+OutFolderName+" --pe1-1 "+OutScriptDict[Locus]['R1_R2'][0]+" --pe1-2 "+OutScriptDict[Locus]['R1_R2'][1]+" -t 4 -m 64\n"
			OutScript2.append(Line)

OutFileName1 = OutFolder+"spades_script_"+OutMode+"1.sh"
OutFileName2 = OutFolder+"spades_script_"+OutMode+"2.sh"

Line = "cat "+OutFileName2+" | parallel --joblog "+OutFolder+"parallel_log.log\n"
OutScript1.append(Line)

OutFile = open(OutFileName1, 'w')
for Line in OutScript1:
	OutFile.write(Line)
OutFile.close()
OutFile = open(OutFileName2, 'w')
for Line in OutScript2:
	OutFile.write(Line)
OutFile.close()


if OutMode == "separate":
	print("The shell script to analyze the individuals separately using spades is %s.\n" % (OutFileName1))
	sys.stderr.write("The shell script to analyze the individuals separately using spades is %s.\n" % (OutFileName1))
elif OutMode == "together":
	print("The shell script to analyze the all individuals together using spades is %s.\n" % (OutFileName1))
	sys.stderr.write("The shell script to analyze the all individuals together using spades is %s.\n" % (OutFileName1))

LocusList = OutScriptDict.keys()
LocusList = sorted(LocusList)

OutFileName = OutFolder+OutFilePre+"Locus_List.txt"
OutFile = open(OutFileName, 'w')
for Locus in LocusList:
	if OutMode == "separate":
		Line = Locus+"\t"+",".join(OutScriptDict[Locus].keys())+"\n"
	elif OutMode == "together":
		Line = Locus+"\ttogether\n"
	OutFile.write(Line)
OutFile.close()

print("The list of loci and the individuals that have sequences for those loci was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The list of loci and the individuals that have sequences for those loci was written to the file %s.\n" % (OutFileName))

