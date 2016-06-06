#! /usr/bin/env python

#trans_bcparse_2reads.py version 1.1 2 March 2015 Abby Moore
#This script is supposed to read the transcriptome files and parse them according to barcode.
#It then makes a separate file for each barcode, labeled according to individual.
#It does not label the sequences according to individual, but leaves them as fastq files.
#I read the whole file, sort it into a dictionary, and then write the output
#to individual files instead of writing the individual files as I go.  But I don't know
#if that's the most efficient way or not??
#Version 1.1 of the script also looks at the barcode on the second read to classify the sequences.
#example:
'''trans_bcparse_2reads.py ~/transcriptomes/TS27_2015_01_26/TS27_1.csv ~/transcriptomes/TS27_2015_01_26/AJM004Bhead.fastq ~/transcriptomes/TS27_2015_01_26/AJM004Bhead_R2.fastq ~/transcriptomes/TS27_2015_01_26/TS27_inds_temp temp_
'''

Usage = '''
trans_bcparse_2reads.py version 1.0
This script reads fastq files and separates them into individuals according to inline barcodes.
trans_bcparse.py [name of file with barcodes] [name of file to parse--read one]
[name of file to parse--read two]
[name of output folder] [prefix for output file]
The barcode file must be tab delimited and have the barcode listed before the individual name.
'''

import sys #We want to be able to get information from the command line
from collections import defaultdict #We want to be able to make dictionaries with multiple levels
import gzip #We want to be able to open gzipped files directly.
from itertools import izip #We want to be able to look at two files simultaneously (without having to
#have them both fully in memory at the same time

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 5:
	print Usage
else:
	BCFileName = sys.argv[1]
	InFileName1 = sys.argv[2]
	InFileName2 = sys.argv[3]
	OutFolder = sys.argv[4]
	OutFilePre = sys.argv[5]

#The major lists and dictionaries we will make:
BCDict = { } #The dictionary of barcodes and their individuals
SeqDict = defaultdict(dict) #The default dictionary of individuals and their sequences
#format: SeqDict[IndName][SeqName=1st line of fastq] = TempSeq=2nd,3rd,4th lines of fastq
SeqDict2 = defaultdict(dict) #The default dictionary of the second reads
#Note that the SeqName is the same is in SeqDict, while the _actual_ first line of the second read
#sequence is included in TempSeq!!
GoodBC = 0 #The number of good barcodes found
BadBC = 0 #The number of bad barcodes found
BCTTT = 0 #The number of sequences whose barcodes consist entirely of Ts
TotalSeqs = 0 #The total number of sequences examined (should be the same as GoodBC + BadBC)

#making sure the output folder ends in a slash
if OutFolder[-1:] != "/":
	OutFolder += "/"

InFile = open(BCFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	BCDict[Line[0]] = Line[1]
InFile.close()

print("%d barcodes were read from the file %s.\n" % (len(BCDict), BCFileName))
sys.stderr.write("%d barcodes were read from the file %s.\n" % (len(BCDict), BCFileName))

BothGood = 0
BothBad = 0
GoodButDifferent = 0
FGood = 0
RGood = 0
if InFileName1.endswith('.gz'):
	InFile1 = gzip.open(InFileName1, 'rU')
else:
	InFile1 = open(InFileName1, 'rU')
if InFileName2.endswith('.gz'):
	InFile2 = gzip.open(InFileName2, 'rU')
else:
	InFile2 = open(InFileName2, 'rU')
LineNum = 0
for Line1, Line2 in izip(InFile1, InFile2):
	Line1 = Line1.strip('\n').strip('\r')
	Line2 = Line2.strip('\n').strip('\r')
	SeqLine = (LineNum + 4) % 4
	if SeqLine == 0: #Code is copied directly (first line of fastq formatted file)
		SeqName = Line1
		TempSeq2 = Line2+"\n"
		LineNum += 1
	elif SeqLine == 1: #Barcode is read and sequence is copied
		Barcode = Line1[0:6]
		TempSeq = Line1[6:]+"\n"
		Barcode2 = Line2[0:6]
		TempSeq2 += Line2[6:]+"\n"
		TotalSeqs += 1
		try: #See if the barcode is in the dictionary
			IndName = BCDict[Barcode]
		except KeyError: #If no valid barcode is found, the individual
		#name is changed to "error".
			IndName = 'error'
			if Barcode == 'TTTTTT':
				BCTTT += 1
		try: #Doing the same for the second barcode
			IndName2 = BCDict[Barcode2]
		except KeyError: #If no valid barcode is found, the individual
		#name is changed to "error".
			IndName2 = 'error'
		if IndName == 'error':
			if IndName2 == 'error':
				BothBad += 1
				BadBC += 1
			else:
				RGood += 1
				IndName = IndName2
				GoodBC += 1
		else:
			if IndName2 == 'error':
				FGood += 1
				GoodBC += 1
			elif IndName2 == IndName2:
				BothGood += 1
				GoodBC += 1
			else:
				GoodButDifferent += 1
				print("ERROR!!!  The sequence has two different barcodes: %s and %s!\n" % (IndName, IndName2))
				sys.stderr.write("ERROR!!!  The sequence has two different barcodes: %s and %s!\n" % (IndName, IndName2))
		LineNum += 1
	elif SeqLine == 2:#The third line of the sequence is not very important, but is copied anyway
		TempSeq = TempSeq + Line1 + "\n"
		TempSeq2 += Line2 + "\n"
		LineNum += 1
	elif SeqLine == 3: #The Phred score is copied without the bases for the barcode
		if IndName != 'error':
			TempSeq += Line1[6:]
			TempSeq2 += Line2[6:]
			SeqDict[IndName][SeqName] = TempSeq
			SeqDict2[IndName][SeqName] = TempSeq2
		LineNum += 1
InFile.close()

#A check to make sure the classifying by barcodes is more or less working:
if (GoodBC+BadBC) != TotalSeqs:
	print("ERROR!!!  %d sequences were read, but only %d were classified according to their barcode.\n" \
	% (TotalSeqs, GoodBC+BadBC))
	sys.stderr.write("ERROR!!!  %d sequences were read, but only %d were classified according to their barcode.\n" \
	% (TotalSeqs, GoodBC+BadBC))

print("%d sequences were read from the file %s.\n" % (TotalSeqs, InFileName1))
sys.stderr.write("%d sequences were read from the file %s.\n" % (TotalSeqs, InFileName1))
print("%d of these had good barcodes and %d of these had bad barcodes.\n" % (GoodBC, BadBC))
sys.stderr.write("%d of these had good barcodes and %d of these had bad barcodes.\n" % (GoodBC, BadBC))
print("%d of the bad barcodes consisted entirely of Ts.\n" % (BCTTT))
sys.stderr.write("%d of the bad barcodes consisted entirely of Ts.\n" % (BCTTT))
print("%d of the sequences with good barcodes had good barcodes on both reads, while %d had good barcodes on one read only.\n" % (BothGood, FGood+RGood))
sys.stderr.write("%d of the sequences with good barcodes had good barcodes on both reads, while %d had good barcodes on one read only.\n" % (BothGood, FGood+RGood))
print("Of the sequences with only one good barcode, %d of these had good barcodes in the forward read only\n" % (FGood))
sys.stderr.write("Of the sequences with only one good barcode, %d of these had good barcodes in the forward read only\n" % (FGood))
print("and %d of these had good barcodes in the rear read only.\n" % (RGood))
sys.stderr.write("and %d of these had good barcodes in the rear read only.\n" % (RGood))

#writing the output files
for IndName in SeqDict.keys():
	IndSeqs = 0
	OutFileName1 = OutFolder+OutFilePre+IndName+"_R1.fastq"
	OutFileName2 = OutFolder+OutFilePre+IndName+"_R2.fastq"
	OutFile1 = open(OutFileName1, 'a')
	OutFile2 = open(OutFileName2, 'a')
	for SeqName in SeqDict[IndName].keys():
		OutLine1 = SeqName+"\n"+SeqDict[IndName][SeqName]+"\n"
		OutFile1.write(OutLine1)
		OutLine2 = SeqDict2[IndName][SeqName] + "\n"
		OutFile2.write(OutLine2)
		IndSeqs += 1
	OutFile1.close()
	OutFile2.close()
	print("%d sequences for the individual %s were written to the files %s and %s.\n" % (IndSeqs, IndName, OutFileName1, OutFileName2))
	sys.stderr.write("%d sequences for the individual %s were written to the files %s and %s.\n" % (IndSeqs, IndName, OutFileName1, OutFileName2))
