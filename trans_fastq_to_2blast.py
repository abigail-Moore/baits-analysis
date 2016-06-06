#! /usr/bin/env python

#trans_fastq_to_2blast.py version 1.1 9 March 2016 Abby Moore
#This script reads gzipped fastq files (ideally those that have been parsed by 
#trans_bcparse.py) and makes them into fasta files 
#for blasting.  It also writes a script to run blastn on those files.
#version 1.0 9 March 2015
#version 1.1 9 March 2016 two modes have been added: parallel for running on a single computer or a single node of a cluster and
#array for running on slurm clusters

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
trans_fastq_to_2blast.py InFolder InFilePre IndListFileName OutFolder OutFilePre BlastDB BlastFilePre Mode[Parallel, Array] NextScript
'''

Usage = '''
trans_fastq_to_2blast.py version 1.1
This script prepares fastq files for blasting.
trans_fastq_to_blast.py [input folder] [input file prefix] [barcode file] 
[folder for output files] [output file prefix] [name of blast database]
[prefix for blast output files]
[mode: either Parallel or Array]
[the next script--if mode is Array]
'''

import sys #We want to be able to get information from the command line
from collections import defaultdict #We want to be able to make dictionaries with multiple levels
import gzip #We want to be able to open gzipped files directly.
import subprocess #We want to be able to talk to the command line.

print ("%s\n" % (" ".join(sys.argv)))

ModeList = ['Parallel', 'Array']

if len(sys.argv) < 9:
	sys.exit("ERROR!!  This script requires 8 or 9 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
InFilePre = sys.argv[2]
IndListFileName = sys.argv[3]
OutFolder = sys.argv[4]
OutFilePre = sys.argv[5]
BlastDB = sys.argv[6]
BlastFilePre = sys.argv[7]
Mode = sys.argv[8]
if (Mode in ModeList) == False:
	sys.exit("ERROR!!  You wanted the mode %s, but it can only be one of the following: %s.\n %s" % (Mode, ", ".join(ModeList), Usage))
if Mode == "Array":
	if len(sys.argv) != 10:
		sys.exit("ERROR!!  When run in Array mode, this script requires 9 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
	NextScript = sys.argv[9]

#############################################################################

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


#############################################################################

IndList = [ ] #List of the individuals for which files will be made
TotalSeqs = 0
GoodInds = 0
BadInds = 0

#making sure the output and input folders end in slashes
if InFolder[-1:] != "/":
	InFolder += "/"
if OutFolder[-1:] != "/":
	OutFolder += "/"

#Getting the names of the individuals from the list of individuals
InFile = open(IndListFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	IndList.append(Line[0])
InFile.close()

OutIndDict = { }#OutIndDict[IndName] = NumIndFiles
for IndName in IndList:
	for Ending in ['_R1','_R2']:
		SeqDict = { }
		#try to open the gzipped file
		NumIndFiles = 1
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
			#counting sequences and files
			IndSeqs = 0
			OutFileName = OutFolder+OutFilePre+IndName+Ending+"_"+str(NumIndFiles)+".fa"
			OutFile = open(OutFileName, 'w')
			for SeqName in SeqDict:
				OutLine = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
				OutFile.write(OutLine)
				IndSeqs += 1
				#we want a new file for every hundred thousand sequences
				if (IndSeqs % 100000) == 0:
					#so close the old file and add its information to the dictionary
					OutFile.close()
					NumIndFiles += 1
					#then start the counter over and open the new file
					OutFileName = OutFolder+OutFilePre+IndName+Ending+"_"+str(NumIndFiles)+".fa"
					OutFile = open(OutFileName, 'w')
					OutLine = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
					OutFile.write(OutLine)
					IndSeqs += 1
			OutFile.close()
			OutIndDict[IndName] = NumIndFiles
			print("%d sequences for the individual %s were written to %d files with names such as %s.\n" % (IndSeqs, IndName, NumIndFiles, OutFileName))
			sys.stderr.write("%d sequences for the individual %s were written to %d files with names such as %s.\n" % (IndSeqs, IndName, NumIndFiles, OutFileName))
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
				OutFileName = OutFolder+OutFilePre+IndName+Ending+"_"+str(NumIndFiles)+".fa"
				OutFile = open(OutFileName, 'w')
				for SeqName in SeqDict:
					OutLine = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
					OutFile.write(OutLine)
					IndSeqs += 1
					#we want a new file for every hundred thousand sequences
					if (IndSeqs % 100000) == 0:
						OutFile.close()
						NumIndFiles += 1
						OutFileName = OutFolder+OutFilePre+IndName+Ending+"_"+str(NumIndFiles)+".fa"
						OutFile = open(OutFileName, 'w')
						OutLine = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
						OutFile.write(OutLine)
						IndSeqs += 1
				OutFile.close()
				OutIndDict[IndName] = NumIndFiles
				print("%d sequences for the individual %s were written to %d files with names such as %s.\n" % (IndSeqs, IndName, NumIndFiles, OutFileName))
				sys.stderr.write("%d sequences for the individual %s were written to %d files with names such as %s.\n" % (IndSeqs, IndName, NumIndFiles, OutFileName))
				GoodInds += 0.5
			except IOError: #If neither the gzipped nor the non-gzipped file exit
				print("There were no sequences for %s.\n" % (IndName))
				sys.stderr.write("There were no sequences for %s.\n" % (IndName))
				BadInds += 0.5

print("In total, %d sequences from %d individuals were read from files with names such as %s.\n" % (TotalSeqs, GoodInds, InFileName))
sys.stderr.write("In total, %d sequences from %d individuals were read from files with names such as %s.\n" % (TotalSeqs, GoodInds, InFileName))
print("%d individuals did not have sequences.\n" % (BadInds))
sys.stderr.write("%d individuals did not have sequences.\n" % (BadInds))

OutFileName = OutFolder+OutFilePre+"Num_Ind_Files_List.txt"
OutFile = open(OutFileName, 'w')
for IndName in sorted(OutIndDict.keys()):
	OutFile.write(IndName+"\t"+str(OutIndDict[IndName])+"\n")
OutFile.close()

if Mode == "Parallel":
	OutScript1 = ["#! /bin/bash\n" ]
	OutScript2 = [ ]
	#each individual's sequences have been broken up into a number of jobs
	for IndName in OutIndDict:
		for Ending in ['_R1_', '_R2_']:
			for FileNum in range(1,OutIndDict[IndName]+1):
				Line = "blastn -db "+BlastDB+" -query "+OutFolder+OutFilePre+IndName+Ending+str(FileNum)+".fa -out "+OutFolder+BlastFilePre+IndName+Ending+str(FileNum)+".out -outfmt '6 std qlen slen' -task blastn\n"
				OutScript2.append(Line)
	OutFileName1 = OutFolder+BlastFilePre+"blast_script1.sh"
	OutFileName2 = OutFolder+BlastFilePre+"blast_script2.sh"
	#and these will be run in parallel
	Line = "cat "+OutFileName2+" | parallel --jobs 4 --joblog "+OutFolder+BlastFilePre+"parallel_log.log\n"
	OutScript1.append(Line)
	OutFileWriting(OutFileName1, OutScript1)
	OutFileWriting(OutFileName2, OutScript2)
	#We do not need to call the script with the rest of the commands in Parallel mode, like we do in Array mode, because it continues on automatically after this.
	print("The script to run blastn in parallel on these fasta files was written to %s.\n" % (OutFileName1))
	sys.stderr.write("The script to run blastn in parallel on these fasta files was written to %s.\n" % (OutFileName1))
	
elif Mode == "Array":
	SBatchList = [ ]
	JobIDPrev = ""
	#each individual will get its own blast array
	for IndName in OutIndDict:
		OutScript = ["#!/bin/bash\n"]
		OutScript.append("#SBATCH -J BLAST_"+IndName.split("_")[-1]+"\n")
		OutScript.append("#SBATCH -t 1:00:00\n#SBATCH --array=1-"+str(OutIndDict[IndName])+"\n")
		OutScript.append("module load blast\nID=$SLURM_ARRAY_TASK_ID\n")
		OutScript.append("blastn -db "+BlastDB+" -query "+OutFolder+OutFilePre+IndName+"_R1_${ID}.fa -out "+OutFolder+BlastFilePre+IndName+"_R1_${ID}.out -outfmt '6 std qlen slen' -task blastn\n")
		OutScript.append("blastn -db "+BlastDB+" -query "+OutFolder+OutFilePre+IndName+"_R2_${ID}.fa -out "+OutFolder+BlastFilePre+IndName+"_R2_${ID}.out -outfmt '6 std qlen slen' -task blastn\n")
		OutFileName =OutFolder+BlastFilePre+IndName+"blscript.sh"
		OutFileWriting(OutFileName, OutScript)
		#if that is the first individual, its array will be run immediately
		if JobIDPrev == "":
			OutLine = "sbatch "+OutFileName
			SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
			JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
		#the array's of subsequent individuals will be run after the previous individual finishes (to avoid using _all_ of the nodes on the server)
		else:
			OutLine = "sbatch -d afterok:"+JobIDPrev+" "+OutFileName
			SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
			JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
		SBatchList.append(JobIDPrev)
	print("%d job arrays were submitted to blast the sequences.\n" % (len(SBatchList)))
	sys.stderr.write("%d job arrays were submitted to blast the sequences.\n" % (len(SBatchList)))
	#this will take place only after the job arrays for all of the individuals have been run.
	OutLine = "sbatch -d afterok:"+":".join(SBatchList)+" "+NextScript
	SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
	JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
	print("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
	sys.stderr.write("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))

