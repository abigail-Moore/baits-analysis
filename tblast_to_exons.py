#! /usr/bin/env python

#tblast_to_exons.py version 1.0 1 April 2015 Abby Moore
#This script reads the output from blasting the spades contigs for a given locus
#against the baits that have been divided into exons to see which exons are present
#in each contig.  Then it chooses the best exons and makes an output file that just 
#has those exons in it.

import sys#to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels

Usage = '''
tblast_to_exons.py version 1.0
The script determines which contigs have which exons and makes a fasta file
containing the chosen contigs.
tblast_to_exons.py [tab-delimited file with blast output (tab) sequences]
[tab-delimited file with gene (tab) # of exons] [directory of blast output files] 
[directory of contig files] [output directory] [output file prefix]
[outfile mode: separate, together, or both][for blasting: folder containing the
sequences][barcode file, with sequence names][sequence file prefix]'''

#example:
'''
tblast_to_exons.py ~/transcriptomes/TS31_1/sandbox/alaAT_file_list.txt ~/transcriptomes/TS31_1/sandbox/alaAT_exon_list.txt ~/transcriptomes/TS31_1/sandbox/ ~/transcriptomes/TS31_1/sandbox/ ~/transcriptomes/TS31_1/sandbox/ b3_ together ~/transcriptomes/TS31_1/sandbox/ ~/transcriptomes/TS31_1/sandbox/bcfiletemp.txt a1_
tblast_to_exons.py BlastListFileName ExonFileName BlastFolder ContigFolder OutFolder OutFilePre OutMode[separate, together, both] SeqFolder BCFileName SeqFilePre
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 11:
	sys.exit("Error! This script requires 10 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
else:
	BlastListFileName = sys.argv[1]
	ExonFileName = sys.argv[2]
	BlastFolder = sys.argv[3]
	ContigFolder = sys.argv[4]
	OutFolder = sys.argv[5]
	OutFilePre = sys.argv[6]
	OutMode = sys.argv[7]
	if (OutMode != 'separate') and (OutMode != 'together') and (OutMode != 'both'):
		sys.exit("Error!  The output mode can only be separate, together, or both, but you wrote %s.\n" % (OutMode))
	SeqFolder = sys.argv[8]
	BCFileName = sys.argv[9]
	SeqFilePre = sys.argv[10]

#adding slashes to the folders, if necessary:
if BlastFolder[-1] != "/":
	BlastFolder += "/"
if ContigFolder[-1] != "/":
	ContigFolder += "/"
if OutFolder[-1] != "/":
	OutFolder += "/"
if SeqFolder[-1] != "/":
	SeqFolder += "/"

#making lists and things to fill out
BlastFileDict = { } #BlastFileDict[BlastFile] = SeqFile
ExonNumDict = { } #ExonNumDict[Locus] = #exons
ContigSeqDict = defaultdict(dict) #ContigSeqDict[Locus][Contig] = ContigSeq
ContigDict = defaultdict(dict) #ContigDict[Locus][Contig] = [list of exons]
ExonContigDict = defaultdict(dict) #ExonContigDict[Locus][Exon] = [list of contigs]
SharedExonDict = defaultdict(dict) #SharedExonDict[Locus][Exon] = [list of other exons that share contigs with that exon]
LocusBlockDict = defaultdict(list) #LocusBlockDict[Locus] = [list of blocks of exons]
ContigsWanted = defaultdict(list) #ContigsWanted[Locus] = [list of contigs]

#reading the list of blast output files
InFile = open(BlastListFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	BlastFileDict[Line[0]] = Line[1]
InFile.close()
print("The names of %d blast output files and their corresponding fasta files were read from the file %s.\n" % (len(BlastFileDict), BlastListFileName))
sys.stderr.write("The names of %d blast output files and their corresponding fasta files were read from the file %s.\n" % (len(BlastFileDict), BlastListFileName))

#reading the list of #s of exons in each gene
InFile = open(ExonFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	ExonNumDict[Line[0]] = int(Line[1])
InFile.close()
print("The numbers of exons for each of %d genes were read from the file %s.\n" % (len(ExonNumDict), ExonFileName))
sys.stderr.write("The numbers of exons for each of %d genes were read from the file %s.\n" % (len(ExonNumDict), ExonFileName))

#reading the blast files and adding that information to the ContigDict:
for BlastFile in BlastFileDict:
	#first, reading the sequences from the corresponding sequence file
	TempSeqDict = { }
	InFileName = ContigFolder+BlastFileDict[BlastFile]
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n')
		if Line[0] == ">":
			Contig = Line[1:]
			TempSeqDict[Contig] = ""
		else:
			TempSeqDict[Contig] += Line
	InFile.close()
	#then reading the blast file
	InFileName = BlastFolder+BlastFile
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		eValue = float(Line[10])
		#if the blast hit is above the threshold e-value
		if eValue < 1e-8:
			Contig = Line[0]
			BlastHit = Line[1].split('+')
			Locus = BlastHit[0]
			Exon = int(BlastHit[1][4:])
			#adding the contig sequence to the dictionary
			ContigSeqDict[Locus][Contig] = TempSeqDict[Contig]
			#recording the fact that that contig blasted to that exon
			try:
				ContigDict[Locus][Contig].append(Exon)
			except KeyError:
				ContigDict[Locus] = defaultdict(list)
				ContigDict[Locus][Contig].append(Exon)
			#recording the fact that that exon had a contig that blasted to it
			try:
				ExonContigDict[Locus][Exon].append(Contig)
			except KeyError:
				ExonContigDict[Locus] = defaultdict(list)
				ExonContigDict[Locus][Exon].append(Contig)
	InFile.close()
#sorting and getting rid of duplicates
for Locus in ContigDict:
	for Contig in ContigDict[Locus]:
		ListTemp = ContigDict[Locus][Contig]
		ContigDict[Locus][Contig] = sorted(list(set(ListTemp)))
	for Exon in ExonContigDict[Locus]:
		ListTemp = ExonContigDict[Locus][Exon]
		ExonContigDict[Locus][Exon] = sorted(list(set(ListTemp)))

#checking to see if all exons are present:
for Locus in ExonNumDict:
	ExonRange = range(1, ExonNumDict[Locus]+1)
	MissingExons = 0
	if ExonContigDict[Locus] == { }:
		print("Locus %s has no contigs.\n" % (Locus))
		sys.stderr.write("Locus %s has no contigs.\n" % (Locus))
	else:
		for Exon in ExonRange:
			if ExonContigDict[Locus][Exon] == [ ]:
				print("Locus %s has no contigs that cover exon %d.\n" % (Locus, Exon))
				#sys.stderr.write("Locus %s has no contigs that cover exon %d.\n" % (Locus, Exon))
				MissingExons += 1
		if MissingExons == 0:
			print("Locus %s has contigs that cover all %d exons.\n" % (Locus, ExonNumDict[Locus]))
			sys.stderr.write("Locus %s has contigs that cover all %d exons.\n" % (Locus, ExonNumDict[Locus]))
		else:
			print("%d of the %d exons of locus %s do not have any contigs.\n" % (MissingExons, ExonNumDict[Locus], Locus))
			sys.stderr.write("%d of the %d exons of locus %s do not have any contigs.\n" % (MissingExons, ExonNumDict[Locus], Locus))
	
#Now I need to determine which blocks the exons fall into.
#first, determining which exons share contigs with which other exons
for Locus in ContigDict:
	SharedExonDict[Locus] = defaultdict(list)
	for Contig in ContigDict[Locus]:
		for Exon in ContigDict[Locus][Contig]:
			SharedExonDict[Locus][Exon] += ContigDict[Locus][Contig]
#then going back through that information and...
for Locus in SharedExonDict:
	TempList = [ ]
	for Exon in SharedExonDict[Locus]:
		#condensing the lists
		ListTemp = SharedExonDict[Locus][Exon]
		ListTemp2 = sorted(list(set(ListTemp)))
		SharedExonDict[Locus][Exon] = sorted(list(set(ListTemp)))
		#perhaps the exon does not share contigs with any other exons
		if len(ListTemp2) == 1:
			LocusBlockDict[Locus].append([Exon])
		#or if it does, determining how large of a block of exons can be made
		elif len(ListTemp2) > 1:
			if (Exon in TempList) == True:
				TempList += ListTemp2
				TempList = list(set(TempList))
			else:
				if TempList != [ ]:
					LocusBlockDict[Locus].append(TempList)
				TempList = ListTemp2
	if TempList != [ ]:
		LocusBlockDict[Locus].append(TempList)


#deciding which contigs we want
for Locus in ContigDict:
	for Contig in ContigDict[Locus]:
		#we want all contigs that cover more than one locus
		#***Do I want this to be able to vary, instead of always being 1?***
		if len(ContigDict[Locus][Contig]) > 1:
			ContigsWanted[Locus].append(Contig)
		else:
			#we also want contigs that cover a single locus, if they are the only contig
			#that includes that locus
			if len(ExonContigDict[Locus][ContigDict[Locus][Contig][0]]) == 1:
				ContigsWanted[Locus].append(Contig)
			#it is also possible that a given locus has multiple contigs, but that they never
			#overlap with those of other loci, because the surrounding intron(s) is/are very long
			elif SharedExonDict[Locus][ContigDict[Locus][Contig][0]] == [ContigDict[Locus][Contig][0]]:
				ContigsWanted[Locus].append(Contig)


#Now writing all of the output files!
LocusList = sorted(ContigDict.keys())
#first, the fasta files of the contigs
OutFileList = [ ]
if (OutMode == 'together') or (OutMode == 'both'):
	OutFileName = OutFolder+OutFilePre+"all.fa"
	OutFileList.append(OutFileName)
	OutFile = open(OutFileName, 'w')
	for Locus in LocusList:
		for Contig in ContigsWanted[Locus]:
			Line = ">"+Locus+"+"+Contig+"\n"+ContigSeqDict[Locus][Contig]+"\n"
			OutFile.write(Line)
	OutFile.close()
	print("All of the desired contig sequences were written to the file %s.\n" % (OutFileName))
	sys.stderr.write("All of the desired contig sequences were written to the file %s.\n" % (OutFileName))
if (OutMode == 'separate') or (OutMode == 'both'):
	for Locus in ContigsWanted:
		OutFileName = OutFolder+OutFilePre+Locus+".fa"
		OutFileList.append(OutFileName)
		OutFile = open(OutFileName, 'w')
		for Contig in ContigsWanted[Locus]:
			Line = ">"+Locus+"+"+Contig+"\n"+ContigSeqDict[Locus][Contig]+"\n"
			OutFile.write(Line)
		OutFile.close()
	print("The contigs from %d loci were written to separate files with names such as %s.\n" % (len(ContigsWanted), OutFileName))
	sys.stderr.write("The contigs from %d loci were written to separate files with names such as %s.\n" % (len(ContigsWanted), OutFileName))

#then the file of which exons each contig includes
OutList = ["Locus\tContig\tExons\tIn_Output\n"]
#ContigDict[Locus][Contig] = [list of exons]
for Locus in LocusList:
	ContigList = sorted(ContigDict[Locus].keys())
	for Contig in ContigList:
		Line = Locus+"\t"+Contig+"\t"+",".join(str(ExonNum) for ExonNum in ContigDict[Locus][Contig])+"\t"
		if (Contig in ContigsWanted[Locus]):
			Line += "yes\n"
		else:
			Line += "no\n"
		OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Exons_per_Contig.txt"
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()
print("The list of exons each locus includes was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The list of exons each locus includes was written to the file %s.\n" % (OutFileName))

#then the files of which exons group together:
OutList = ["Locus\t#Blocks\tExons_per_Block\n"]
for Locus in LocusList:
	Line = Locus+"\t"+str(len(LocusBlockDict[Locus]))
	for LocusBlock in LocusBlockDict[Locus]:
		Line += "\t"+(",").join(str(LocusNum) for LocusNum in LocusBlock)
	Line += "\n"
	OutList.append(Line)
OutFileName = OutFolder+OutFilePre+"Exon_Blocks.txt"
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()
print("Information on the blocks of exons into which the contigs for each locus were divided was written to the file %s.\n" % (OutFileName))
sys.stderr.write("Information on the blocks of exons into which the contigs for each locus were divided was written to the file %s.\n" % (OutFileName))

#then making the script to write the blast databases
OutList = ['#! /bin/bash\n']
BlastFileList = [ ]
for OutFileName in OutFileList:
	Line = 'makeblastdb -in '+OutFileName+' -out '+OutFileName[:-3]+' -dbtype nucl\n'
	OutList.append(Line)
	BlastFileList.append(OutFileName[:-3]+"\n")
OutFileName = OutFolder+OutFilePre+'BlastDBscript.sh'
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()
print("The script to make blast databases from these fasta files is %s.\n" % (OutFileName))
sys.stderr.write("The script to make blast databases from these fasta files is %s.\n" % (OutFileName))

#then making a script to run blast.
InFile = open(BCFileName, 'rU')
IndList = [ ]
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	IndList.append(Line[1])
if (OutMode == 'together') or (OutMode == 'both'):
	OutList = ['#! /bin/bash\n']
	for Ind in IndList:
		for Ending in ['_R1', '_R2']:
			Line = 'blastn -db '+OutFolder+OutFilePre+"all -query "+SeqFolder+SeqFilePre+Ind+Ending+".fa -out "+OutFolder+OutFilePre+Ind+Ending+".out -outfmt '6 std qlen slen' -task blastn\n"
			OutList.append(Line)
OutFileName = OutFolder+OutFilePre+'Blastn_script.sh'
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()
print("The script to run blast was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The script to run blast was written to the file %s.\n" % (OutFileName))
