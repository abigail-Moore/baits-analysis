#! /usr/bin/env python

#tassembly_to_blast.py version 1.0 3 June 2015 Abby Moore
#This script parses the output from running spades on the sequences that
#blasted to individual loci (using a shell script written by tblast_to_spades.py --together)
#to make a blast database from these loci and to blast the original sequences back against them.
#The LocusFile (prefix_Locus_List.txt written by tblast_to_spades.py) should have two columns:
#Locus [0]
#"together" instead of a list of individuals [1]

import sys#We want to be able to talk to the command line
import subprocess#We want to be able to get input from shell processes
from collections import defaultdict #I am sure that these will come in handy somewhere
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences

#Example:
'''
tassembly_to_blast.py InFolder LocusFileName SeqFolder SeqFilePre NumIndSeqsFileName BlastDBName OutFolder OutFilePre Mode[Parallel, Array] NextScript
'''

print ("%s\n" % (" ".join(sys.argv)))

Usage ='''
tassembly_to_loci.py version 1.0
This is a script to combine the masurca contigs from multiple individuals for a
single locus so they can be analyzed further.
tassembly_to_loci.py
[folder with the folders of sequences arranged according to locus]
[file with the list of loci] 
[Folder with the original sequences] 
[prefix for the sequence files] 
[tab-delimitted file produced by trans_fastq_to_2blast.py that has the
IndividualName [tab] Number of Sequence Files]
[name of fasta file for the blast database to which the contig sequences should 
be added]
[output folder] 
[output file prefix, which is the same for the new blast database and the files 
made from blasting the sequences against that database]
[mode: either Parallel or Array]
[script for the next analysis, if in Array mode]
'''

ModeList = ['Parallel', 'Array']

if len(sys.argv) < 10:
	sys.exit("ERROR!!  This script requires 9 or 10 additional arguments, and you provided %d.  %s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
LocusFileName = sys.argv[2]
SeqFolder = sys.argv[3]
SeqFilePre = sys.argv[4]
if SeqFilePre == "none":
	SeqFilePre = ""
NumIndSeqsFileName = sys.argv[5]
BlastDBName = sys.argv[6]
OutFolder = sys.argv[7]
OutFilePre = sys.argv[8]
if OutFilePre == "none":
	OutFilePre = ""
Mode = sys.argv[9]
if (Mode in ModeList) == False:
	sys.exit("ERROR!!  You wanted the mode %s, but it can only be one of the following: %s.\n %s" % (Mode, ", ".join(ModeList), Usage))
if Mode == "Array":
	if len(sys.argv) != 11:
		sys.exit("ERROR!!  When run in Array mode, this script requires 10 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
	NextScript = sys.argv[10]

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


LocusDict = { } #List of loci and their associated individuals, read from LocusFileName
ContigDict = defaultdict(dict)#Dictionary of contigs arranged according to locus
#and individual: ContigDict[Locus][Ind][SeqName] = Seq
ContigNumDict = { }#The number of contigs per locus
OutIndDict = DictFromFile(NumIndSeqsFileName, 0, 1)#OutIndDict[IndName] = NumFiles
IndList = sorted(OutIndDict.keys())
print("%d individuals will be examined.\n" % (len(IndList)))
sys.stderr.write("%d individuals will be examined.\n" % (len(IndList)))

NumSeqs = 0

if InFolder[-1:] != "/":
	InFolder += "/"
if OutFolder[-1:] != "/":
	OutFolder += "/"
if SeqFolder[-1:] != "/":
	SeqFolder += "/"

#Getting the names of the loci and the individuals that have those loci from the locus file
InFile = open(LocusFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	Locus = Line[0]
	LocusDict[Locus] = Line[1].split(",")
InFile.close()
print("The list of %d loci was read from the file %s.\n" % (len(LocusDict), LocusFileName))
sys.stderr.write("The list of %d loci was read from the file %s.\n" % (len(LocusDict), LocusFileName))

#getting the sequences for that locus from the various contigs.fasta files
IndNumDict = defaultdict(int)
for Locus in LocusDict:
	ContigDict[Locus] = defaultdict(dict)
	NumSeqs = 0
	for Ind in LocusDict[Locus]:
		InFileExists = "no"
		#see if the a file exists for that individual for that locus for that format:
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
		if InFileExists == "yes":
			for record in SeqIO.parse(InFileName, 'fasta'):
				#only take sequences of a certain length, but I don't know what this should be*****
				if len(str(record.seq)) > 200:
					SeqName = Locus+"+"+ContigPre+"+"+record.id
					ContigDict[Locus][Ind][SeqName] = str(record.seq)
					NumSeqs += 1
					IndNumDict[Ind] += 1
	ContigNumDict[Locus] = NumSeqs
	print("%d contigs were found in a total of %d individuals for locus %s.\n" % (NumSeqs, len(ContigDict[Locus].keys()), Locus))
	sys.stderr.write("%d contigs were found in a total of %d individuals for locus %s.\n" % (NumSeqs, len(ContigDict[Locus].keys()), Locus))

#writing everything to the outfile
OutFileName = OutFolder+OutFilePre+"post_spades_combined.fa"
OutList = [ ]
for Locus in ContigDict.keys():
	if len(ContigDict[Locus].keys()) > 0:
		NumSeqs = 0
		for Ind in ContigDict[Locus].keys():
			for SeqName in ContigDict[Locus][Ind].keys():
				NumSeqs += 1
				OutLine = ">"+SeqName+"\n"+ContigDict[Locus][Ind][SeqName]+"\n"
				OutList.append(OutLine)
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

#writing the blast script:
if Mode == "Parallel":
	BlastList1 = ['#! /bin/bash\n']
	BlastList2 = [ ]
	Line = "cat "+BlastDBName+".fa "+OutFileName+" > "+OutFolder+OutFilePre+"spades_baits_combined.fa\n"
	Line += "makeblastdb -in "+OutFolder+OutFilePre+"spades_baits_combined.fa -out "+OutFolder+OutFilePre+"spades_baits_combined -dbtype nucl\n"
	BlastList1.append(Line)
	for IndName in OutIndDict:
		for Ending in ['_R1_', '_R2_']:
			NumFiles = int(OutIndDict[IndName])
			for FileNum in range(1,NumFiles+1):
				Line = "blastn -db "+OutFolder+OutFilePre+"spades_baits_combined -query "+SeqFolder+SeqFilePre+IndName+Ending+str(FileNum)+".fa -out "+OutFolder+OutFilePre+IndName+Ending+str(FileNum)+".out -outfmt '6 std qlen slen' -task blastn\n"
				BlastList2.append(Line)
	#Now I need to write the script for blasting the other files against this database.
	OutFileName1 = OutFolder+OutFilePre+"BlastScript1.sh"
	OutFileName2 = OutFolder+OutFilePre+"BlastScript2.sh"
	#running blast in parallel
	Line = "cat "+OutFileName2+" | parallel --joblog "+OutFolder+OutFilePre+"parallel_log.log\n"
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
	OutCode = subprocess.call("cat "+BlastDBName+".fa "+OutFolder+OutFilePre+"post_spades_combined.fa > "+OutFolder+OutFilePre+"spades_baits_combined.fa", shell = True, stdout=subprocess.PIPE)
	OutCode = subprocess.call("makeblastdb -in "+OutFolder+OutFilePre+"spades_baits_combined.fa -out "+OutFolder+OutFilePre+"spades_baits_combined -dbtype nucl", shell = True, stdout=subprocess.PIPE)
	print("The blast database %s has hopefully been made.  If it has not, you will get a lot of scripts that will never start running.....\n" % (OutFolder+OutFilePre+"spades_baits_combined"))
	sys.stderr.write("The blast database %s has hopefully been made.  If it has not, you will get a lot of scripts that will never start running.....\n" % (OutFolder+OutFilePre+"spades_baits_combined"))
	JobIDPrev = ""
	SBatchList = [ ]
	for IndName in OutIndDict:
		OutScript = ["#!/bin/bash\n"]
		OutScript.append("#SBATCH -J BLAST_"+IndName.split("_")[-1]+"\n")
		OutScript.append("#SBATCH -t 1:30:00\n#SBATCH --array=1-"+str(OutIndDict[IndName])+"\n")
		OutScript.append("module load blast\nID=$SLURM_ARRAY_TASK_ID\n")
		OutScript.append("blastn -db "+OutFolder+OutFilePre+"spades_baits_combined -query "+SeqFolder+SeqFilePre+IndName+"_R1_${ID}.fa -out "+OutFolder+OutFilePre+IndName+"_R1_${ID}.out -outfmt '6 std qlen slen' -task blastn\n")
		OutScript.append("blastn -db "+OutFolder+OutFilePre+"spades_baits_combined -query "+SeqFolder+SeqFilePre+IndName+"_R2_${ID}.fa -out "+OutFolder+OutFilePre+IndName+"_R2_${ID}.out -outfmt '6 std qlen slen' -task blastn\n")
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
