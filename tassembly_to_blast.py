#! /usr/bin/env python

#tassembly_to_blast.py version 1.0 3 June 2015 Abby Moore
#This script parses the output from running spades on the sequences that
#blasted to individual loci (using a shell script written by tblast_to_spades.py --together)
#to make a blast database from these loci and to blast the original sequences back against them.
#The LocusFile (prefix_Locus_List.txt written by tblast_to_spades.py) should have two columns:
#Locus [0]
#"together" instead of a list of individuals [1]

import sys#We want to be able to talk to the command line
from collections import defaultdict #I am sure that these will come in handy somewhere
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences

#Example:
'''
tassembly_to_blast.py InFolder LocusFileName SeqFolder SeqFilePre IndListFileName BlastDBName OutFolder OutFilePre
'''

print ("%s\n" % (" ".join(sys.argv)))

Usage ='''
tassembly_to_loci.py version 1.0
This is a script to combine the masurca contigs from multiple individuals for a
single locus so they can be analyzed further.
tmasurca_to_loci.py [folder with the folders of sequences arranged according to locus]
[file with the list of loci] [Folder with the original sequences] [prefix for
the sequence files] [barcode file--list of individuals in 2nd column] [name of
fasta file for the blast database to which the contig sequences should be added]
[output folder] [output file prefix, which is the same for the new blast 
database and the files made from blasting the sequences against that database]'''

if len(sys.argv) != 9:
	sys.exit(Usage)
else:
	InFolder = sys.argv[1]
	LocusFileName = sys.argv[2]
	SeqFolder = sys.argv[3]
	SeqFilePre = sys.argv[4]
	if SeqFilePre == "none":
		SeqFilePre = ""
	IndListFileName = sys.argv[5]
	BlastDBName = sys.argv[6]
	OutFolder = sys.argv[7]
	if sys.argv[8] == "none":
		OutFilePre = ""
	else:
		OutFilePre = sys.argv[8]

LocusDict = { } #List of loci and their associated individuals, read from LocusFileName
ContigDict = defaultdict(dict)#Dictionary of contigs arranged according to locus
#and individual: ContigDict[Locus][Ind][SeqName] = Seq
ContigNumDict = { }#The number of contigs per locus
IndList = [ ] #List of individuals to be read from IndListFileName
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

#getting the list of individuals to be blasted against these sequences
InFile = open(IndListFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	IndList.append(Line[0])
InFile.close()
print("The list of %d individuals was read from the file %s.\n" % (len(IndList), IndListFileName))
sys.stderr.write("The list of %d individuals was read from the file %s.\n" % (len(IndList), IndListFileName))

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

NumSeqsDict = defaultdict(list)
for Ind in IndNumDict:
	NumSeqsDict[IndNumDict[Ind]].append(Ind)
IndListOut = [ ]
for NumSeqs in sorted(NumSeqsDict.keys(), reverse = True):
	IndListOut += NumSeqsDict[NumSeqs]

#writing the blast script:
BlastList1 = ['#! /bin/bash\n']
BlastList2 = [ ]
Line = "cat "+BlastDBName+" "+OutFileName+" > "+OutFolder+OutFilePre+"spades_baits_combined.fa\n"
Line += "makeblastdb -in "+OutFolder+OutFilePre+"spades_baits_combined.fa -out "+OutFolder+OutFilePre+"spades_baits_combined -dbtype nucl\n"
BlastList1.append(Line)
for IndName in IndListOut:
	Line = "blastn -task blastn -db "+OutFolder+OutFilePre+"spades_baits_combined -query "+SeqFolder+SeqFilePre+IndName+"_R1.fa -out "+OutFolder+OutFilePre+IndName+"_R1.out -outfmt '6 std qlen slen'\n"
	Line += "blastn -task blastn -db "+OutFolder+OutFilePre+"spades_baits_combined -query "+SeqFolder+SeqFilePre+IndName+"_R2.fa -out "+OutFolder+OutFilePre+IndName+"_R2.out -outfmt '6 std qlen slen'\n"
	BlastList2.append(Line)
#Now I need to write the script for blasting the other files against this database.
OutFileName1 = OutFolder+OutFilePre+"BlastScript1.sh"
OutFileName2 = OutFolder+OutFilePre+"BlastScript2.sh"
Line = "cat "+OutFileName2+" | parallel --jobs = "+NCores+" --joblog "+OutFolder+OutFilePre+"parallel_log.log\n"
OutFileName1.append(Line)

OutFile = open(OutFileName1,'w')
for Line in BlastList1:
	OutFile.write(Line)
OutFile.close()
OutFile = open(OutFileName2,'w')
for Line in BlastList2:
	OutFile.write(Line)
OutFile.close()
print("The script for making a blast database out of these sequences is %s.\n" % (OutFileName1))
sys.stderr.write("The script for making a blast database out of these sequences is %s.\n" % (OutFileName1))

