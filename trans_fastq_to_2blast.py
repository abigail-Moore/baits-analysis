#! /usr/bin/env python

#trans_fastq_to_2blast.py version 1.0 9 March 2015 Abby Moore
#This script reads gzipped fastq files (ideally those that have been parsed by 
#trans_bcparse.py) and makes them into fasta files 
#for blasting.  It also writes a script to run blastn on those files.

#Decoding the first line of fastq format
# @HWI-ST558:60:C00B3ACXX:7:1206:8335:28258 1:N:0:ATCACG
# It has the following parts (with meaning according to Wikipedia):
# @HWI-ST558: the unique instrument name (same for all runs): 0
# 60: the run id (same within each run): 1
# C00B3ACXX: the flowcell id (same within each run): 2
# 7: flowcell lane (same within each sample): 3
# 1206: tile number within the flowcell lane (changes within a sample) **need this**: 4
# 8335: 'x'-coordinate of the cluster within the tile (changes within a sample) **need this**: 5
# 28258: 'y'-coordinate of the cluster within the tile (changes within a sample) **need this**: 6:0
# 1: the member of a pair, 1 or 2 (paired-end or mate-pair reads only): 6:1
# N: Y if the read fails filter (read is bad), N otherwise (I assume we do not have any Y samples.): 7
# 0: 0 when none of the control bits are on, otherwise it is an even number: 8
# ATCACG: index sequence: 9 [This isn't always present, depending on whether or not a third-read barcode was used.]

#example:
'''
trans_fastq_to_2blast.py ~/transcriptomes/TS27_2015_01_26/TS27_2015_01_26_inds/ t1_ ~/transcriptomes/TS27_2015_01_26/TS27_1.csv ~/transcriptomes/TS27_2015_01_26/TS27_2015_01_26_inds/ a1_
trans_fastq_to_2blast.py InFolder InFilePre BCFileName OutFolder OutFilePre BlastDB BlastFilePre
'''

Usage = '''
trans_fastq_to_2blast.py version 1.0
This script prepares fastq files for blasting.
trans_fastq_to_blast.py [input folder] [input file prefix] [barcode file] 
[folder for output files] [output file prefix] [name of blast database]
[prefix for blast output files]
'''

import sys #We want to be able to get information from the command line
from collections import defaultdict #We want to be able to make dictionaries with multiple levels
import gzip #We want to be able to open gzipped files directly.

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 8:
	sys.exit("ERROR!!  This script requires 7 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
InFilePre = sys.argv[2]
BCFileName = sys.argv[3]
OutFolder = sys.argv[4]
OutFilePre = sys.argv[5]
BlastDB = sys.argv[6]
BlastFilePre = sys.argv[7]

IndList = [ ] #List of the individuals for which files will be made
TotalSeqs = 0
GoodInds = 0
BadInds = 0
OutFileList = [ ] #List of files that the program writes

#making sure the output and input folders end in slashes
if InFolder[-1:] != "/":
	InFolder += "/"
if OutFolder[-1:] != "/":
	OutFolder += "/"

#Getting the names of the individuals from the barcodes file
InFile = open(BCFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	IndList.append(Line[1])
InFile.close()

for IndName in IndList:
	for Ending in ['_R1','_R2']:
		SeqDict = { }
		try:
			InFileName = InFolder + InFilePre + IndName + Ending + ".fastq.gz"
			InFile = gzip.open(InFileName, 'rU')
			LineNum = 0
			for Line in InFile:
				Line = Line.strip('\n').strip('\r')
				SeqLine = (LineNum + 4) % 4
				if SeqLine == 0: #Code is copied directly (first line of fastq formatted file)
					Line = Line.split(" ")[0].split(":")
					SeqName = IndName+"_"+Line[4]+"_"+Line[5]+"_"+Line[6]+Ending
				elif SeqLine == 1: #Sequence is added to the dictionary
					Seq = Line
					SeqDict[SeqName] = Seq
					TotalSeqs += 1
				LineNum += 1		
			InFile.close()		
			IndSeqs = 0
			OutFileName = OutFolder+OutFilePre+IndName+Ending+".fa"
			OutFile = open(OutFileName, 'w')
			for SeqName in SeqDict:
				OutLine = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
				OutFile.write(OutLine)
				IndSeqs += 1
			OutFile.close()
			OutFileList.append(OutFileName)
			print("%d sequences for the individual %s were written to the file %s.\n" % (IndSeqs, IndName, OutFileName))
			sys.stderr.write("%d sequences for the individual %s were written to the file %s.\n" % (IndSeqs, IndName, OutFileName))
			GoodInds += 0.5
		except IOError:
			try: #Trying to open the non-gzipped file
				InFileName = InFolder + InFilePre + IndName + Ending + ".fastq"
				InFile = open(InFileName, 'rU')
				LineNum = 0
				for Line in InFile:
					Line = Line.strip('\n').strip('\r')
					SeqLine = (LineNum + 4) % 4
					if SeqLine == 0: #Code is copied directly (first line of fastq formatted file)
						Line = Line.split(" ")[0].split(":")
						SeqName = IndName+"_"+Line[4]+"_"+Line[5]+"_"+Line[6]+Ending
					elif SeqLine == 1: #Sequence is added to the dictionary
						Seq = Line
						SeqDict[SeqName] = Seq
						TotalSeqs += 1
					LineNum += 1		
				InFile.close()		
				IndSeqs = 0
				OutFileName = OutFolder+OutFilePre+IndName+Ending+".fa"
				OutFile = open(OutFileName, 'w')
				for SeqName in SeqDict:
					OutLine = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
					OutFile.write(OutLine)
					IndSeqs += 1
				OutFile.close()
				OutFileList.append(OutFileName)
				print("%d sequences for the individual %s were written to the file %s.\n" % (IndSeqs, IndName, OutFileName))
				sys.stderr.write("%d sequences for the individual %s were written to the file %s.\n" % (IndSeqs, IndName, OutFileName))
				GoodInds += 0.5
			except IOError: #If neither the gzipped nor the non-gzipped file exit
				print("There were no sequences for %s.\n" % (IndName))
				sys.stderr.write("There were no sequences for %s.\n" % (IndName))
				BadInds += 0.5

print("In total, %d sequences from %d individuals were read from files with names such as %s.\n" % (TotalSeqs, GoodInds, InFileName))
sys.stderr.write("In total, %d sequences from %d individuals were read from files with names such as %s.\n" % (TotalSeqs, GoodInds, InFileName))
print("%d individuals did not have sequences.\n" % (BadInds))
sys.stderr.write("%d individuals did not have sequences.\n" % (BadInds))

BlastFileList = [ ]
FastaFileList = [ ]
Line = "#! /bin/bash\n"
BlastFileList.append(Line)
for IndName in IndList:
	for Ending in ['_R1', '_R2']:
		Line = "blastn -db "+BlastDB+" -query "+OutFolder+OutFilePre+IndName+Ending+".fa -out "+OutFolder+BlastFilePre+IndName+Ending+".out -outfmt '6 std qlen slen' -task blastn\n"
		BlastFileList.append(Line)
		FastaFileList.append(OutFolder+OutFilePre+IndName+Ending+".fa\n")
OutFileName = OutFolder+BlastFilePre+"blast_script.sh"
OutFile = open(OutFileName, 'w')
for Line in BlastFileList:
	OutFile.write(Line)
OutFile.close()
print("The script to run blastn on these fasta files was written to %s.\n" % (OutFileName))
sys.stderr.write("The script to run blastn on these fasta files was written to %s.\n" % (OutFileName))

#***This part hasn't been specifically tested***
OutFileName = OutFolder+OutFilePre+"fasta_file_list.txt"
OutFile = open(OutFileName, 'w')
for Line in FastaFileList:
	OutFile.write(Line)
OutFile.close()
print("The list of fasta files was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of fasta files was written to %s.\n" % (OutFileName))
