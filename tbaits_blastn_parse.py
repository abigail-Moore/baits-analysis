#! /usr/bin/env python

#tbaits_blast_parse.py version 1.1 3 June 2015 Abby Moore
#This script reads the blast output from blasting the captured sequences (using
#blastn) against the blast database made from the bait sequences.
#It requires a tab-delimitted list of locus names.  These can either be the same 
#name in both columns or the name that is in the blast output in the first column
#and a modified name in the second column (e.g., asp2	asp)
#It writes various files, among others, a list of sequences with blast hits to the
#various loci ([OutFilePre]Seqs_to_Loci.txt), which will be read by tblast_to_fastq.py
#this version has been modified so that it no longer looks for plastid and nrDNA sequences

import sys #to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels

#Example:
'''
tbaits_blastn_parse.py InFolderBl NumIndSeqsFileName LocusFileName BLFilePre OutFolder OutFilePre
'''

Usage = '''
tbaits_blastn_parse.py version 1.1
This script analyzes the output from blasting the captured sequences of each individual
against the blast database made from the bait sequences.
tbaits_blast_parse.py
[folder containing blast output]
[tab-delimitted file produced by trans_fastq_to_2blast.py that has the
IndividualName [tab] Number of Sequence Files]
[file containing the updated locus names]
[prefix for blast output files]
[output folder] 
[prefix for output files]
'''
print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 7:
	sys.exit("Error! tbaits_blastn_parse.py requires six additional arguments and you supplied %d.  %s" % (len(sys.argv)-1, Usage))
InFolderBl = sys.argv[1]
NumIndSeqsFileName = sys.argv[2]
LocusFileName = sys.argv[3]
BLFilePre = sys.argv[4]
OutFolder = sys.argv[5]
OutFilePre = sys.argv[6]

#################################################################################

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

###################################################################################
#making sure the output and input folders end in slashes
if InFolderBl[-1:] != "/":
	InFolderBl += "/"
if OutFolder[-1:] != "/":
	OutFolder += "/"

OutIndDict = DictFromFile(NumIndSeqsFileName, 0, 1)#OutIndDict[IndName] = NumFiles
IndList = sorted(OutIndDict.keys())
print("%d individuals will be examined.\n" % (len(IndList)))
sys.stderr.write("%d individuals will be examined.\n" % (len(IndList)))
LocusDict = { } #The dictionary of the locus names from the blast files and the final output files for those names
#LocusDict[LocusName] = Locus
SeqBlastDict = defaultdict(dict) #The dictionary of sequence names and their blast hits
BlastHitsDict = defaultdict(dict) #The dictionary of the blast hits, so we know which individuals had which blast hits.
BlastHitsCounts = defaultdict(dict) #The dictionary of the counts of blast hits per individual per locus
IndHitsDict = defaultdict(dict) #The dictionary of the number of blast hits per individual
IndSumDict = defaultdict(dict) #The dictionary with the summary of the blast hits
BlastIDDict = defaultdict(int)#The dictionary that shows how often each combination of blast hits was found.
MaxHits = 0
EndingList = ['_R1_','_R2_']

#Getting the names of the loci from the locus file
InFile = open(LocusFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	LocusDict[Line[0]] = Line[1]
InFile.close()

#looking at the blast output
SeqBlastDict = defaultdict(dict)
for IndName in IndList:
	for Ending in EndingList:
		for FileNum in range(1,int(OutIndDict[IndName])+1):
			InFileName = InFolderBl+BLFilePre+IndName+Ending+str(FileNum)+".out"
			InFile = open(InFileName, 'rU')
			try:
				for Line in InFile:
					Line = Line.strip('\n').strip('\r').split('\t')
					SeqName = Line[0][:-3]
					BlHit = Line[1].split("+")[0]
					Locus = LocusDict[BlHit]
					eValue = float(Line[10])
					#*****Here is where you change the eValue cutoff*****
					if eValue < 1e-16:
						try:
							SeqBlastDict[IndName][SeqName][Ending].append(Locus)
						except KeyError:
							SeqBlastDict[IndName][SeqName] = defaultdict(list)
							SeqBlastDict[IndName][SeqName][Ending].append(Locus)
			except SystemError:
				sys.stderr.write("SystemError!!")
			InFile.close()

print("Blast results were read from %d output files with names of the form %s.\n" % (len(IndList)*2, InFileName))
sys.stderr.write("Blast results were read from %d output files with names of the form %s.\n" % (len(IndList)*2, InFileName))

#condensing the lists
for IndName in IndList:
	#first, set the summary dictionary to zero:
	IndSumDict[IndName]['Bait'] = 0
	#Then, going through the sequences for that individual
	for SeqName in SeqBlastDict[IndName].keys():
		#removing multiple references to the same blast hit
		for Ending in EndingList:
			SeqBlastDict[IndName][SeqName][Ending] = list(set(SeqBlastDict[IndName][SeqName][Ending]))
			SeqBlastDict[IndName][SeqName]['Total'] += SeqBlastDict[IndName][SeqName][Ending]
		SeqBlastDict[IndName][SeqName]['Total'] = list(set(SeqBlastDict[IndName][SeqName]['Total']))
		ListTemp = SeqBlastDict[IndName][SeqName]['Total']
		#adding information about the number of hits that individual had to the dictionary
		NumHits = len(ListTemp)
		if NumHits > MaxHits:
			MaxHits = NumHits
		try:
			IndHitsDict[IndName][NumHits] += 1
		except KeyError:
			IndHitsDict[IndName][NumHits] = 1
		#calculating more things, if that sequence had blast hits
		if NumHits != 0:
			IndSumDict[IndName]['Bait'] += 1
			for BlHit in SeqBlastDict[IndName][SeqName]['Total']:
				try:
					BlastHitsDict[BlHit][IndName][NumHits].append(SeqName)
				except KeyError:
					BlastHitsDict[BlHit][IndName] = defaultdict(list)
					BlastHitsDict[BlHit][IndName][NumHits] = [SeqName]
				if NumHits == 1:
					TempName = IndName+"_1"
				elif NumHits > 1:
					TempName = IndName+"_m"
				try:
					BlastHitsCounts[BlHit][TempName] += 1
				except KeyError:
					BlastHitsCounts[BlHit][TempName] = 1
			SeqBlastID = "_".join(sorted(SeqBlastDict[IndName][SeqName]['Total']))
			BlastIDDict[SeqBlastID] += 1

#Writing the results to files that make some kind of sense:
BlHitList = sorted(BlastHitsCounts.keys())
OutList = [ ]
Head = "LocusName"
for IndName in IndList:
	Head += "\t"+IndName+"_1\t"+IndName+"_m"
OutList.append(Head)
for BlHit in BlHitList:
	Line = "\n"+BlHit
	for IndName in IndList:
		try:
			Line += "\t"+str(BlastHitsCounts[BlHit][IndName+"_1"])
		except KeyError:
			Line += "\t"
		try:
			Line += "\t"+str(BlastHitsCounts[BlHit][IndName+"_m"])
		except KeyError:
			Line += "\t"
	OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Hits_per_Locus.txt"
OutFileWriting(OutFileName, OutList)

print("Information about the number of hits for each locus from each individual was written to %s.\n" % (OutFileName))
sys.stderr.write("Information about the number of hits for each locus from each individual was written to %s.\n" % (OutFileName))

OutList = [ ]
#SeqStatsList = [ ]
#StatsHead = "Name\tTotal_Reads\tReads_with_no_hits\tReads_with_hits\n"
#SeqStatsList.append(StatsHead)
HitRange = range(0,MaxHits+1,1)
Head = "Hits_per_Sequence"
for Num in HitRange:
	Head += "\t"+str(Num)
OutList.append(Head)
for IndName in IndList:
	SeqswithHits = 0
	Line = "\n"+IndName
	for Num2 in HitRange:
		try:
			Line += "\t"+str(IndHitsDict[IndName][Num2])
			if (IndHitsDict[IndName][Num2] > 0) and (Num2 > 0):
				SeqswithHits += IndHitsDict[IndName][Num2]
		except KeyError:
			Line += "\t0"
	OutList.append(Line)
	#StatsLine = IndName+"\t"+str(SeqswithHits+IndHitsDict[IndName][0])+"\t"+str(IndHitsDict[IndName][0])+"\t"+str(SeqswithHits)+"\t"
	#SeqStatsList.append(StatsLine)
OutFileName = OutFolder+OutFilePre+"Seqs_with_Hits.txt"
OutFileWriting(OutFileName, OutList)

print("Information about the number of sequences each individual had that had a certain number of hits was written to the file %s.\n" % (OutFileName))
sys.stderr.write("Information about the number of sequences each individual had that had a certain number of hits was written to the file %s.\n" % (OutFileName))

#This is taken out because, since we no longer read the sequence file, we no longer know how many sequences there were.
'''
OutFileName = OutFolder+OutFilePre+"Hit_Distribution.txt"
OutFile = open(OutFileName, 'w')
for Line in SeqStatsList:
	OutFile.write(Line)
OutFile.close()

print("Information about the distribution of these hits across the genome was written to the file %s.\n" % (OutFileName))
sys.stderr.write("Information about the distribution of these hits across the genome was written to the file %s.\n" % (OutFileName))
'''

#BlastHitsDict[BlHit][IndName][NumHits] = [SeqName]
OutList = [ ]
Head = "Locus_Name\tIndividual_Name\tNumber_of_Hits\tSequence_Names"
OutList.append(Head)
for BlHit in BlHitList:
	for IndName in BlastHitsDict[BlHit].keys():
		for NumHits in BlastHitsDict[BlHit][IndName].keys():
			Line = "\n" + BlHit + "\t" + IndName + "\t" + str(NumHits) + "\t" + ";".join(BlastHitsDict[BlHit][IndName][NumHits])
			OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Seqs_to_Loci.txt"
OutFileWriting(OutFileName, OutList)

print("The list of sequences that go with each locus was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of sequences that go with each locus was written to %s.\n" % (OutFileName))

#BlastIDDict[BlastID] = NumHits
OutList = ["BlastID\tNumber_of_Hits\n"]
for BlastID in sorted(BlastIDDict.keys()):
	Line = BlastID+"\t"+str(BlastIDDict[BlastID])+"\n"
	OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Hit_Combinations.txt"
OutFileWriting(OutFileName, OutList)

print("The list of combinations of blast hits and the number of times each combination was found was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of combinations of blast hits and the number of times each combination was found was written to %s.\n" % (OutFileName))
