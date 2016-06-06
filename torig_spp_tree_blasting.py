#! /usr/bin/env python

#torig_spp_tree_blasting.py version 1.0 14 April 2016 Abby Moore
#This script blasts files of sequences for different individuals against a 
#pre-made blast database.
#The sequence files can be numbered IndName_FilePost_1-... from a list
#modified from tassembly_to_blast.py

import sys#We want to be able to talk to the command line
import subprocess#We want to be able to get input from shell processes

#Example:
'''
torig_spp_tree_blasting.py SeqFolder SeqFilePre NumIndSeqsFileName BlastDBName OutFolder OutFilePre Mode[Parallel, Array] NextScript
'''

print ("%s\n" % (" ".join(sys.argv)))

Usage ='''
torig_spp_tree_blasting.py version 1.0
This is a script to blast files against a pre-made blast database.
[Folder with the original sequences] 
[prefix for the sequence files] 
[tab-delimitted file produced by trans_fastq_to_2blast.py that has the
IndividualName [tab] Number of Sequence Files]
[name of the blast database]
[output folder] 
[output file prefix, which is the same for the new blast database and the files 
made from blasting the sequences against that database]
[mode: either Parallel or Array]
[script for the next analysis, if in Array mode]
'''

ModeList = ['Parallel', 'Array']

if len(sys.argv) < 8:
	sys.exit("ERROR!!  This script requires 7 or 8 additional arguments, and you provided %d.  %s" % (len(sys.argv)-1, Usage))
SeqFolder = sys.argv[1]
SeqFilePre = sys.argv[2]
if SeqFilePre == "none":
	SeqFilePre = ""
NumIndSeqsFileName = sys.argv[3]
BlastDBName = sys.argv[4]
OutFolder = sys.argv[5]
OutFilePre = sys.argv[6]
if OutFilePre == "none":
	OutFilePre = ""
Mode = sys.argv[7]
if (Mode in ModeList) == False:
	sys.exit("ERROR!!  You wanted the mode %s, but it can only be one of the following: %s.\n %s" % (Mode, ", ".join(ModeList), Usage))
if Mode == "Array":
	if len(sys.argv) != 9:
		sys.exit("ERROR!!  When run in Array mode, this script requires 10 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
	NextScript = sys.argv[8]

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

#DictFromFile makes a dictionary from a tab-delimited file, where the keys are in columns KC, and
#the values are in column VC
#from tcontig_selection.py, which was a modified version of the function in from tbaits_intron_removal.py
def DictFromFile(FileName, KC, VC):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[KC]] = Line[VC]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is OutIndDict
	
#############################################################################


OutIndDict = DictFromFile(NumIndSeqsFileName, 0, 1)#OutIndDict[IndName] = NumFiles
IndList = sorted(OutIndDict.keys())
print("%d individuals will be examined.\n" % (len(IndList)))
sys.stderr.write("%d individuals will be examined.\n" % (len(IndList)))

if OutFolder[-1:] != "/":
	OutFolder += "/"
if SeqFolder[-1:] != "/":
	SeqFolder += "/"

#writing the blast script:
if Mode == "Parallel":
	BlastList1 = ['#! /bin/bash\n']
	BlastList2 = [ ]
	for IndName in OutIndDict:
		for Ending in ['_R1_', '_R2_']:
			for FileNum in range(1,int(OutIndDict[IndName])+1):
				Line = "blastn -db "+BlastDBName+" -query "+SeqFolder+SeqFilePre+IndName+Ending+str(FileNum)+".fa -out "+OutFolder+OutFilePre+IndName+Ending+str(FileNum)+".out -outfmt '6 std qlen slen' -task blastn\n"
				BlastList2.append(Line)
	#Now I need to write the script for blasting the other files against this database.
	OutFileName1 = OutFolder+OutFilePre+"BlastScript1.sh"
	OutFileName2 = OutFolder+OutFilePre+"BlastScript2.sh"
	#running blast in parallel
	Line = "cat "+OutFileName2+" | parallel --jobs 4 --joblog "+OutFolder+OutFilePre+"parallel_log.log\n"
	BlastList1.append(Line)
	OutFile = open(OutFileName1,'w')
	for Line in BlastList1:
		OutFile.write(Line)
	OutFile.close()
	OutFile = open(OutFileName2,'w')
	for Line in BlastList2:
		OutFile.write(Line)
	OutFile.close()
	print("The script for making a blast database out of these sequences and blasting the old sequences against in parallel is %s.\n" % (OutFileName1))
	sys.stderr.write("The script for making a blast database out of these sequences and blasting the old sequences against in parallel is %s.\n" % (OutFileName1))
elif Mode == "Array":
	JobIDPrev = ""
	SBatchList = [ ]
	for IndName in OutIndDict:
		OutScript = ["#!/bin/bash\n"]
		OutScript.append("#SBATCH -J BLAST_"+IndName.split("_")[-1]+"\n")
		OutScript.append("#SBATCH -t 1:30:00\n#SBATCH --array=1-"+str(OutIndDict[IndName])+"\n")
		OutScript.append("module load blast\nID=$SLURM_ARRAY_TASK_ID\n")
		OutScript.append("blastn -db "+BlastDBName+" -query "+SeqFolder+SeqFilePre+IndName+"_R1_${ID}.fa -out "+OutFolder+OutFilePre+IndName+"_R1_${ID}.out -outfmt '6 std qlen slen' -task blastn\n")
		OutScript.append("blastn -db "+BlastDBName+" -query "+SeqFolder+SeqFilePre+IndName+"_R2_${ID}.fa -out "+OutFolder+OutFilePre+IndName+"_R2_${ID}.out -outfmt '6 std qlen slen' -task blastn\n")
		OutFileName =OutFolder+OutFilePre+IndName+"blscript.sh"
		OutFileWriting(OutFileName, OutScript)
		if JobIDPrev == "":
			OutLine = "sbatch "+OutFileName
			print OutLine
			SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
			JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
		else:
			OutLine = "sbatch -d afterok:"+JobIDPrev+" "+OutFileName
			print OutLine
			SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
			JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
		SBatchList.append(JobIDPrev)
	print("%d job arrays were submitted to blast the sequences.\n" % (len(SBatchList)))
	sys.stderr.write("%d job arrays were submitted to blast the sequences.\n" % (len(SBatchList)))
	OutLine = "sbatch -d afterok:"+":".join(SBatchList)+" "+NextScript
	print OutLine
	SBatchOut = subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
	JobIDPrev = SBatchOut.strip('\r').strip('\n').split(" ")[-1]
	print("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
	sys.stderr.write("The subsequent script will be run with the job id %s, when the previous scripts have finished.\n" % (JobIDPrev))
