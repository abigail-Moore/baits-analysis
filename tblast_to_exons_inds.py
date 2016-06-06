#! /usr/bin/env python

#tblast_to_exons.py version 1.0 1 April 2015 Abby Moore
#This script reads the output from blasting the spades contigs for a given locus
#against the baits that have been divided into exons to see which exons are present
#in each contig.  Then it chooses the best exons and makes an output file that just 
#has those exons in it.

import sys#to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels
from Bio.Seq import Seq #to edit sequences
from Bio.Alphabet import IUPAC

Usage = '''
tblast_to_exons.py version 1.0
The script determines which contigs have which exons and makes a fasta file
containing the chosen contigs.
tblast_to_exons.py [tab-delimited file with blast output (tab) sequences]
[tab-delimited file with gene (tab) # of exons] [directory of blast output files] 
[directory of contig files] [output directory] [output file prefix]
[for phylogenetic analysis: folder containing the template alignments (with genomic
sequence on the last line and its cds on the penultimate line][template file
prefix or "none"] [template file ending or "none"] [sequences mode: either
all sequences (all) or just the best sequences (best)] [output mode: either whole
sequences (contig) or exons only (exon)] [exon mode: none (if whole contigs are 
output), separate (separate fasta files for each exon), bycontig (one fasta sequence per
contig, with introns removed), byind (one fasta sequence per individual, with exons
from different contigs combined), byindlist (multiple fasta sequences per individual,
for multiple alleles, as determined by the list) or bylist (separate fasta files 
for larger chunks of each locus)] [exon list file name, if exon mode is bylist]'''

#example:
'''
tblast_to_exons_inds.py ~/transcriptomes/TS31_1/sandbox/alaAT_file_list2.txt ~/transcriptomes/TS31_1/sandbox/alaAT_exon_list.txt ~/transcriptomes/TS31_1/sandbox/ ~/transcriptomes/TS31_1/sandbox/ ~/transcriptomes/TS31_1/sandbox/ s5_ ~/transcriptomes/baits_Bv/ none _Bv_al all exon byind 
tblast_to_exons_inds.py BlastListFileName ExonFileName BlastFolder ContigFolder OutFolder OutFilePre TAlFolder TAlPre TAlPost SeqsMode [all, best] OutMode [contig, exon] ExonMode [none, separate, bycontig, byind, byindlist, bylist] ExonListFile (if ExonMode == bylist)
'''

OutModeList = ['contig', 'exon']
ExonModeList = ['none', 'separate', 'bycontig', 'byind', 'byindlist', 'bylist']
SeqsModeList = ['all', 'best']

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) < 13:
	sys.exit("Error! This script requires 11-12 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
else:
	BlastListFileName = sys.argv[1]
	ExonFileName = sys.argv[2]
	BlastFolder = sys.argv[3]
	ContigFolder = sys.argv[4]
	OutFolder = sys.argv[5]
	OutFilePre = sys.argv[6]
	TAlFolder = sys.argv[7]
	TAlPre = sys.argv[8]
	if TAlPre == "none":
		TAlPre = ""
	TAlPost = sys.argv[9]
	if TAlPost == "none":
		TAlPost = ""
	SeqsMode = sys.argv[10]
	if (SeqsMode in SeqsModeList) == False:
		sys.exit("ERROR!  The sequence mode must be one of the following: %s, but you wrote %s.\n" % (" ".join(SeqsModeList), SeqsMode))
	OutMode = sys.argv[11]
	if (OutMode in OutModeList) == False:
		sys.exit("ERROR!  The output mode must be one of the following: %s, but you wrote %s.\n" % (" ".join(OutModeList), OutMode))
	ExonMode = sys.argv[12]
	if (ExonMode in ExonModeList) == False:
		sys.exit("ERROR!  The exon mode must be one of the following: %s, but you wrote %s.\n" % (" ".join(ExonModeList), ExonMode))
	if (OutMode == 'together') and (ExonMode != 'none'):
		print("ERROR!  The exon mode %s is incompatible with the output mode 'together', so the exon mode will be ignored.\n" % (ExonMode))
		sys.stderr.write("ERROR!  The exon mode %s is incompatible with the output mode 'together', so the exon mode will be ignored.\n" % (ExonMode))
	if ExonMode == "byind":
		print("WARNING!  Exons from different contigs from the same individual will be combined.  You should only use this mode if you are certain that the sequences are single copy!!!\n")
		sys.stderr.write("WARNING!  Exons from different contigs from the same individual will be combined.  You should only use this mode if you are certain that the sequences are single copy!!!\n")
	if (ExonMode == "bylist") or (ExonMode == "byindlist"):
		try:
			ExonListFile = sys.argv[13]
		except IndexError:
			sys.exit("ERROR!  If the exon mode is bylist or byindlist, you also need to supply a file showing how the contigs should be divided.")

#adding slashes to the folders, if necessary:
if BlastFolder[-1] != "/":
	BlastFolder += "/"
if ContigFolder[-1] != "/":
	ContigFolder += "/"
if OutFolder[-1] != "/":
	OutFolder += "/"
if TAlFolder[-1] != "/":
	TAlFolder += "/"

#making lists and things to fill out
BlastFileDict = { } #BlastFileDict[BlastFile] = SeqFile
ExonNumDict = { } #ExonNumDict[Locus] = #exons
ContigSeqDict = defaultdict(dict) #ContigSeqDict[Locus][Ind][Contig] = ContigSeq
ContigDict = defaultdict(dict) #ContigDict[Locus][Ind][Contig] = [list of exons]
ExonContigDict = defaultdict(dict) #ExonContigDict[Locus][Ind][Exon] = [list of contigs]
ExonPosDict = defaultdict(dict) #ExonPosDict[Locus][Ind][Exon][variousinformation] = variousinformation
SharedExonDict = defaultdict(dict) #SharedExonDict[Locus][Ind][Exon] = [list of other exons that share contigs with that exon]
LocusBlockDict = defaultdict(list) #LocusBlockDict[Locus] = [list of blocks of exons]
ContigsWanted = defaultdict(dict) #ContigsWanted[Locus][Ind] = [list of contigs]
OutExonDict = defaultdict(dict) #OutExonDict[Locus][Exon] = [fasta-formatted exons in list format] (only used with some modes)
ExonsWanted = defaultdict(dict) #ExonsWanted[Locus][Ind][Exon] = Contig (only used in some modes)

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

#reading the file with the instructions on how to parse the sequences
if (OutMode == 'exon') and (ExonMode == 'byindlist'): #ContigGroupDict[Locus][Ind][Group][Exon]=list of contigs
	ContigGroupDict = defaultdict(dict)
	InFile = open(ExonListFile,'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		Locus = Line[0]
		Contig = Line[1]
		Ind = Line[2]
		Exon = int(Line[4])
		Group = Line[5]
		try:
			ContigGroupDict[Locus][Ind][Group][Exon].append(Contig)
		except KeyError:
			try:
				ContigGroupDict[Locus][Ind][Group] = defaultdict(list)
				ContigGroupDict[Locus][Ind][Group][Exon].append(Contig)
			except KeyError:
				ContigGroupDict[Locus][Ind] = defaultdict(dict)
				ContigGroupDict[Locus][Ind][Group] = defaultdict(list)
				ContigGroupDict[Locus][Ind][Group][Exon].append(Contig)
	

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
			Ind = Contig.split('+')[0]
			if (int(Line[8]) > int(Line[9])):
				Direc = 'rev'
			else:
				Direc = 'normal'
			#adding the full contig sequence to the dictionary
			try:
				ContigSeqDict[Locus][Ind][Contig]['Whole'] = TempSeqDict[Contig]
			except KeyError:
				ContigSeqDict[Locus][Ind] = defaultdict(dict)
				ContigSeqDict[Locus][Ind][Contig]['Whole'] = TempSeqDict[Contig]
			TempStart = int(Line[6])
			TempEnd = int(Line[7])
			#recording the fact that that contig blasted to that exon
			try:
				ContigDict[Locus][Ind][Contig].append(Exon)
			except KeyError:
				ContigDict[Locus][Ind] = defaultdict(list)
				ContigDict[Locus][Ind][Contig].append(Exon)
			#recording the fact that that exon had a contig that blasted to it
			try:
				ExonContigDict[Locus][Ind][Exon].append(Contig)
			except KeyError:
				ExonContigDict[Locus][Ind] = defaultdict(list)
				ExonContigDict[Locus][Ind][Exon].append(Contig)
			#recording the start and end of the exon
			try:
				if ExonPosDict[Locus][Contig][Exon]['Start'] > TempStart:
					ExonPosDict[Locus][Contig][Exon]['Start'] = TempStart
				if ExonPosDict[Locus][Contig][Exon]['End'] < TempEnd:
					ExonPosDict[Locus][Contig][Exon]['End'] = TempEnd
			except KeyError:
				try:
					ExonPosDict[Locus][Ind][Contig][Exon]['Start'] = TempStart
					ExonPosDict[Locus][Ind][Contig][Exon]['End'] = TempEnd
					ExonPosDict[Locus][Ind][Contig][Exon]['Direc'] = Direc
				except KeyError:
					try:
						ExonPosDict[Locus][Ind][Contig][Exon] = defaultdict(dict)
						ExonPosDict[Locus][Ind][Contig][Exon]['Start'] = TempStart
						ExonPosDict[Locus][Ind][Contig][Exon]['End'] = TempEnd
						ExonPosDict[Locus][Ind][Contig][Exon]['Direc'] = Direc
					except KeyError:
						ExonPosDict[Locus][Ind] = defaultdict(dict)
						ExonPosDict[Locus][Ind][Contig][Exon] = defaultdict(dict)
						ExonPosDict[Locus][Ind][Contig][Exon]['Start'] = TempStart
						ExonPosDict[Locus][Ind][Contig][Exon]['End'] = TempEnd
						ExonPosDict[Locus][Ind][Contig][Exon]['Direc'] = Direc
	InFile.close()

#sorting and getting rid of duplicates
for Locus in ContigDict:
	for Ind in ContigDict[Locus]:
		for Contig in ContigDict[Locus][Ind]:
			ListTemp = ContigDict[Locus][Ind][Contig]
			ContigDict[Locus][Ind][Contig] = sorted(list(set(ListTemp)))
		for Exon in ExonContigDict[Locus][Ind]:
			ListTemp = ExonContigDict[Locus][Ind][Exon]
			ExonContigDict[Locus][Ind][Exon] = sorted(list(set(ListTemp)))
#adding the exon sequences to the ContigSeqDict
for Locus in ExonPosDict:
	for Ind in ExonPosDict[Locus]:
		for Contig in ExonPosDict[Locus][Ind]:
			TempSeq = ContigSeqDict[Locus][Ind][Contig]['Whole']
			for Exon in ExonPosDict[Locus][Ind][Contig]:
				ExonSeq = TempSeq[(ExonPosDict[Locus][Ind][Contig][Exon]['Start']-1):(ExonPosDict[Locus][Ind][Contig][Exon]['End']-1)]
				if ExonPosDict[Locus][Ind][Contig][Exon]['Direc'] == 'rev':
					ExonSeq2 = Seq(ExonSeq, IUPAC.unambiguous_dna)
					ExonSeq = str(ExonSeq2.reverse_complement())
				ContigSeqDict[Locus][Ind][Contig][Exon] = ExonSeq
			if ExonPosDict[Locus][Ind][Contig][Exon]['Direc'] == 'rev':
				TempSeq2 = Seq(TempSeq, IUPAC.unambiguous_dna)
				TempSeq = str(TempSeq2.reverse_complement())
				ContigSeqDict[Locus][Ind][Contig]['Whole'] = TempSeq
				

#Now I need to determine which blocks the exons fall into.
#first, determining which exons share contigs with which other exons
for Locus in ContigDict:
	for Ind in ContigDict[Locus]:
		SharedExonDict[Locus][Ind] = defaultdict(list)
		for Contig in ContigDict[Locus][Ind]:
			for Exon in ContigDict[Locus][Ind][Contig]:
				SharedExonDict[Locus][Ind][Exon] += ContigDict[Locus][Ind][Contig]
#then going back through that information and...
for Locus in SharedExonDict:
	for Ind in SharedExonDict[Locus]:
		TempList = [ ]
		for Exon in SharedExonDict[Locus][Ind]:
			#condensing the lists
			ListTemp = SharedExonDict[Locus][Ind][Exon]
			ListTemp2 = sorted(list(set(ListTemp)))
			SharedExonDict[Locus][Ind][Exon] = sorted(list(set(ListTemp)))

#deciding which contigs we want
for Locus in ContigDict:
	ContigsWanted[Locus] = defaultdict(list)
	for Ind in ContigDict[Locus]:
		for Contig in ContigDict[Locus][Ind]:
			#we want all contigs that cover more than one locus
			#***Do I want this to be able to vary, instead of always being 1?***
			if SeqsMode == 'best':
				if len(ContigDict[Locus][Ind][Contig]) > 1:
					ContigsWanted[Locus][Ind].append(Contig)
				else:
					#we also want contigs that cover a single locus, if they are the only contig
					#that includes that locus
					if len(ExonContigDict[Locus][Ind][ContigDict[Locus][Ind][Contig][0]]) == 1:
						ContigsWanted[Locus][Ind].append(Contig)
					#it is also possible that a given locus has multiple contigs, but that they never
					#overlap with those of other loci, because the surrounding intron(s) is/are very long
					elif SharedExonDict[Locus][Ind][ContigDict[Locus][Ind][Contig][0]] == [ContigDict[Locus][Ind][Contig][0]]:
						ContigsWanted[Locus][Ind].append(Contig)
			elif SeqsMode == 'all':
				ContigsWanted[Locus][Ind].append(Contig)

#Now writing all of the output files!
LocusList = sorted(ContigDict.keys())

#first, the actual sequence files:
#if the whole contigs are output (both exons and introns)
if (OutMode == 'contig'):
	for Locus in ContigsWanted:
		OutFileName = OutFolder+OutFilePre+Locus+".fa"
		OutFile = open(OutFileName, 'w')
		for Ind in ContigsWanted[Locus]:
			for Contig in ContigsWanted[Locus][Ind]:
				Line = ">"+Contig+"_"+Locus+"\n"+ContigSeqDict[Locus][Ind][Contig]['Whole']+"\n"
				OutFile.write(Line)
		OutFile.close()
	print("The contigs from %d loci were written to separate files with names such as %s.\n" % (len(ContigsWanted), OutFileName))
	sys.stderr.write("The contigs from %d loci were written to separate files with names such as %s.\n" % (len(ContigsWanted), OutFileName))
#OutExonDict = defaultdict(dict) #OutExonDict[Locus][Exon] = [fasta-formatted exons in list format]
#separate fasta files for each exon:
elif (OutMode == 'exon') and (ExonMode == 'separate'):
	ExonLocusList = [ ]
	for Locus in ContigsWanted:
		OutExonDict[Locus] = defaultdict(list)
		for Ind in ContigsWanted[Locus]:
			for Contig in ContigsWanted[Locus][Ind]:
				ContigNT = Contig.split("+")
				SeqName = "N"+ContigNT[2].split("_")[1]+"_"+ContigNT[0]+"_"+ContigNT[1]+"+"+"_".join(str(ExonNum) for ExonNum in ExonPosDict[Locus][Ind][Contig].keys())
				for Exon in ContigSeqDict[Locus][Ind][Contig]:
					if Exon != 'Whole':
						TempFasta = ">"+str(Exon)+SeqName+"\n"+ContigSeqDict[Locus][Ind][Contig][Exon]+"\n"
						OutExonDict[Locus][Exon].append(TempFasta)
						Line = Locus+"\t"+Contig+"\t"+Ind+"\t"+SeqName+"\t"+str(Exon)+"\n"
						ExonLocusList.append(Line)
	for Locus in OutExonDict:
		for Exon in OutExonDict[Locus]:
			OutFileName = OutFolder+OutFilePre+str(Exon)+"_"+Locus+".fa"
			OutFile = open(OutFileName, 'w')
			for TempFasta in OutExonDict[Locus][Exon]:
				OutFile.write(TempFasta)
			OutFile.close()
	print("Separate fasta files for each exon from each locus were written, with names such as %s.\n" % (OutFileName))
	sys.stderr.write("Separate fasta files for exon exon from each locus were written, with names such as %s.\n" % (OutFileName))
	OutFileName = OutFolder+OutFilePre+"ExonLocusList.txt"
	OutFile = open(OutFileName, 'w')
	for Line in ExonLocusList:
		OutFile.write(Line)
	OutFile.close()
	print("The list of the contig and sequence names was written to %s.\n" % (OutFileName))
	sys.stderr.write("The list of the contig and sequence names was written to %s.\n" % (OutFileName))
#one fasta file per contig that includes only the exons:
elif (OutMode == 'exon') and (ExonMode == 'bycontig'):
	for Locus in ContigsWanted:
		OutContigs = [ ]
		for Ind in ContigsWanted[Locus]:
			for Contig in ContigsWanted[Locus][Ind]:
				TempExonList = [ ]
				ContigNT = Contig.split("+")
				Line = ">"+"N"+ContigNT[2].split("_")[1]+"+"+"_".join(str(ExonNum) for ExonNum in ExonPosDict[Locus][Ind][Contig].keys())+"+"+ContigNT[0]+"+"+ContigNT[1]+"\n"
				for Exon in ContigSeqDict[Locus][Ind][Contig]:
					if Exon != 'Whole':
						TempExonList.append(str(ContigSeqDict[Locus][Ind][Contig][Exon]))
				Line += "--".join(TempExonList)+"\n"
				OutContigs.append(Line)
		OutFileName = OutFolder+OutFilePre+Locus+"_exons.fa"
		OutFile = open(OutFileName, 'w')
		for Line in OutContigs:
			OutFile.write(Line)
		OutFile.close()
	print("The exons from the contigs of these loci were written to a separate file for each locus with names such as %s.\n" % (OutFileName))
	sys.stderr.write("The exons from the contigs of these loci were written to a separate file for each locus with names such as %s.\n" % (OutFileName))
#one sequence per individual with exons from different contigs
#ExonsWanted = defaultdict(dict) #ExonsWanted[Locus][Ind][Exon] = Contig
#ExonContigDict = defaultdict(dict) #ExonContigDict[Locus][Ind][Exon] = [list of contigs]
elif (OutMode == 'exon') and (ExonMode == 'byind'):
	ExonSourceDict = defaultdict(dict) #ExonSourceDict[Locus][Ind][Exon] = Contig
	for Locus in ExonContigDict:
		OutContigs = [ ]
		ExonSourceDict[Locus] = defaultdict(dict)
		for Ind in ExonContigDict[Locus]:
			TempExonSeqs = [ ]
			TempExonList = sorted(ExonContigDict[Locus][Ind].keys())
			for Exon in TempExonList:
				TempContig = ""
				TempLength = 0
				for Contig in ExonContigDict[Locus][Ind][Exon]:
					if len(ContigSeqDict[Locus][Ind][Contig][Exon]) > TempLength:
						TempContig = Contig
						TempSeq = str(ContigSeqDict[Locus][Ind][Contig][Exon])
						TempLength = len(ContigSeqDict[Locus][Ind][Contig][Exon])
				TempExonSeqs.append(TempSeq)
				ExonSourceDict[Locus][Ind][Exon] = TempContig
			Line = ">"+Ind+"_"+Locus+"+"+"_".join(str(Exon) for Exon in TempExonList)+"\n"+"--".join(TempExonSeqs)+"\n"
			OutContigs.append(Line)
		OutFileName = OutFolder+OutFilePre+Locus+"_byind.fa"
		OutFile = open(OutFileName, 'w')
		for Line in OutContigs:
			OutFile.write(Line)
		OutFile.close()
	print("One file per locus, containing a single sequence per individual, using the best sequence for each exon taken from multiple contigs, was written, with a name such as %s.\n" % (OutFileName))
	sys.stderr.write("One file per locus, containing a single sequence per individual, using the best sequence for each exon taken from multiple contigs, was written, with a name such as %s.\n" % (OutFileName))
elif (OutMode == 'exon') and (ExonMode == 'byindlist'):
	#ContigGroupDict[Locus][Ind][Group][Exon]=list of contigs
	ExonSourceDict = defaultdict(dict) #ExonSourceDict[Locus][Ind][Group][Exon] = Contig
	for Locus in ExonContigDict:
		OutContigs = [ ]
		try:
			for Ind in ContigGroupDict[Locus]:
				ExonSourceDict[Locus][Ind] = defaultdict(dict)
				for Group in ContigGroupDict[Locus][Ind]:
					if Group != "?":
						TempExonList = sorted(ContigGroupDict[Locus][Ind][Group].keys())
						TempExonSeqs = [ ]
						for Exon in TempExonList:
							TempContig = ""
							TempLength = 0
							for Contig in ContigGroupDict[Locus][Ind][Group][Exon]:
								if len(ContigSeqDict[Locus][Ind][Contig][Exon]) > TempLength:
									TempContig = Contig
									TempSeq = str(ContigSeqDict[Locus][Ind][Contig][Exon])
									TempLength = len(ContigSeqDict[Locus][Ind][Contig][Exon])
							TempExonSeqs.append(TempSeq)
							ExonSourceDict[Locus][Ind][Group][Exon] = TempContig
						Line = ">"+Ind+"_"+Locus+"_"+str(Group)+"+"+"_".join(str(Exon) for Exon in TempExonList)+"\n"+"--".join(TempExonSeqs)+"\n"
						OutContigs.append(Line)
		except KeyError:		
			for Ind in ExonContigDict[Locus]:
				TempExonSeqs = [ ]
				TempExonList = sorted(ExonContigDict[Locus][Ind].keys())
				for Exon in TempExonList:
					TempContig = ""
					TempLength = 0
					for Contig in ExonContigDict[Locus][Ind][Exon]:
						if len(ContigSeqDict[Locus][Ind][Contig][Exon]) > TempLength:
							TempContig = Contig
							TempSeq = str(ContigSeqDict[Locus][Ind][Contig][Exon])
							TempLength = len(ContigSeqDict[Locus][Ind][Contig][Exon])
					TempExonSeqs.append(TempSeq)
					ExonSourceDict[Locus][Ind][Exon] = TempContig
				Line = ">"+Ind+"_"+Locus+"+"+"_".join(str(Exon) for Exon in TempExonList)+"\n"+"--".join(TempExonSeqs)+"\n"
				OutContigs.append(Line)
		OutFileName = OutFolder+OutFilePre+Locus+"_byindlist.fa"
		OutFile = open(OutFileName, 'w')
		for Line in OutContigs:
			OutFile.write(Line)
		OutFile.close()
	print("One file per locus, containing one or paralogous sequences per individual, using the best sequence for each exon taken from multiple contigs, was written, with a name such as %s.\n" % (OutFileName))
	sys.stderr.write("One file per locus, containing one or paralogous sequences per individual, using the best sequence for each exon taken from multiple contigs, was written, with a name such as %s.\n" % (OutFileName))
elif (OutMode == 'exon') and (ExonMode == 'bylist'):
	sys.exit("Sorry, we lied, this option is not, in fact, available.")
	
#then the file of which exons each contig includes
OutList = ["Locus\tContig\tExons\tIn_Output\n"]
#ContigDict[Locus][Ind][Contig] = [list of exons]
for Locus in LocusList:
	IndList = sorted(ContigDict[Locus].keys())
	for Ind in IndList:
		ContigList = sorted(ContigDict[Locus][Ind].keys())
		for Contig in ContigList:
			Line = Locus+"\t"+Contig+"\t"+",".join(str(ExonNum) for ExonNum in ContigDict[Locus][Ind][Contig])+"\t"
			if (Contig in ContigsWanted[Locus][Ind]):
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

#if necessary, the file showing which exon came from which contig
#ExonSourceDict[Locus][Ind][Exon] = Contig
if (OutMode == 'exon') and (ExonMode == 'byind'):
	OutList = ["Locus\tIndividual\tExon\tContig_Name\n"]
	for Locus in ExonSourceDict:
		for Ind in ExonSourceDict[Locus]:
			for Exon in ExonSourceDict[Locus][Ind]:
				Line = Locus+"\t"+Ind+"\t"+str(Exon)+"\t"+ExonSourceDict[Locus][Ind][Exon]+"\n"
				OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+"Exon_Source.txt"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	print("The file showing which contig was the source of which exon was written to %s.\n" % (OutFileName))
	sys.stderr.write("The file showing which contig was the source of which exon was written to %s.\n" % (OutFileName))
elif (OutMode == 'exon') and (ExonMode == 'byindlist'):
	OutList = ["Locus\tGroup\tIndividual\tExon\tContig_Name\n"]
	for Locus in ExonSourceDict:
		for Ind in ExonSourceDict[Locus]:
			for Group in ExonSourceDict[Locus][Ind]:
				for Exon in ExonSourceDict[Locus][Ind][Group]:
					Line = Locus+"\t"+Group+"\t"+Ind+"\t"+str(Exon)+"\t"+ExonSourceDict[Locus][Ind][Group][Exon]+"\n"
					OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+"Exon_Source.txt"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	print("The file showing which contig was the source of which exon was written to %s.\n" % (OutFileName))
	sys.stderr.write("The file showing which contig was the source of which exon was written to %s.\n" % (OutFileName))


#then making the script to do further analyses
OutList = ['#! /bin/bash\n']
if (OutMode == 'exon') and (ExonMode == 'separate'):
	for Locus in OutExonDict:
		for Exon in OutExonDict[Locus]:
			Line = "rm "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs.fa\n"
			Line += "rm "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs_al.fa\n"
			Line += "rm "+OutFolder+"RAxML*"+OutFilePre+str(Exon)+"_"+Locus+"\n"
			Line += "cat "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+".fa "+TAlFolder+TAlPre+str(Exon)+"_"+Locus+TAlPost+".fa > "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs.fa\n"
			Line += "mafft --localpair --maxiterate 1000 "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs.fa > "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs_al.fa --quiet\n"
			Line += "./fasta_to_phylip.py "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs_al.fa\n"
			Line += "raxmlHPC -s "+OutFolder+OutFilePre+str(Exon)+"_"+Locus+"_allseqs_al.phy -n "+OutFilePre+str(Exon)+"_"+Locus+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -w "+OutFolder+"\n"
			OutList.append(Line)	
	OutFileName = OutFolder+OutFilePre+"Exon_Alignment.sh"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	print("The script for aligning the exon contigs to the original alignments is %s.\n" % (OutFileName))
	sys.stderr.write("The script for aligning the exon contigs to the original alignments is %s.\n" % (OutFileName))
	print("WARNING!!!  This script will overwrite existing alignments and trees with these prefixes!!\n")
	sys.stderr.write("WARNING!!!  This script will overwrite existing alignments and trees with these prefixes!!\n")
elif (OutMode == 'exon') and (ExonMode == 'bylist'):
	"do something else that is yet to be determined"
else:
	if (OutMode == 'contig'):
		OutFilePost = '_allseqs'
	elif (ExonMode == 'bycontig'):
		OutFilePost = '_exons'
	elif (ExonMode == 'byind'):
		OutFilePost = '_byind'
	elif (ExonMode == 'byindlist'):
		OutFilePost = '_byindlist'
	for Locus in LocusList:
		Line = "rm "+OutFolder+OutFilePre+Locus+OutFilePost+"_all.fa\n"
		Line += "rm "+OutFolder+OutFilePre+Locus+OutFilePost+"_all_al.fa\n"
		Line += "rm "+OutFolder+"RAxML*"+OutFilePre+Locus+OutFilePost+"\n"
		Line += "cat "+OutFolder+OutFilePre+Locus+OutFilePost+".fa "+TAlFolder+TAlPre+Locus+TAlPost+".fa > "+OutFolder+OutFilePre+Locus+OutFilePost+"_all.fa\n"
		Line += "mafft --localpair --maxiterate 1000 "+OutFolder+OutFilePre+Locus+OutFilePost+"_all.fa > "+OutFolder+OutFilePre+Locus+OutFilePost+"_all_al.fa --quiet\n"
		Line += "./fasta_to_phylip.py "+OutFolder+OutFilePre+Locus+OutFilePost+"_all_al.fa\n"
		Line += "raxmlHPC -s "+OutFolder+OutFilePre+Locus+OutFilePost+"_all_al.phy -n "+OutFilePre+Locus+OutFilePost+" -m GTRCAT -p 1234 -f a -N 100 -x 1234 -w "+OutFolder+"\n"
		OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+"Contig_Alignment.sh"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	print("The script for aligning the contigs to the original alignments is %s.\n" % (OutFileName))
	sys.stderr.write("The script for aligning the contigs to the original alignments is %s.\n" % (OutFileName))
	print("WARNING!!!  This script will overwrite existing alignments and trees with these prefixes!!\n")
	sys.stderr.write("WARNING!!!  This script will overwrite existing alignments and trees with these prefixes!!\n")
