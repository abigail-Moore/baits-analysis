#! /usr/bin/env python

#tblast_to_fasta_qual.py version 1.0 27 March 2015
#This script reads the Seqs_to_Loci.txt output 
#from tbaits_blastn_parse.py or tbaits_blast_parse.py
#and writes files for each locus for each individual.
#It expects file with the header:
#Locus_Name	Individual_Name	Number_of_Hits	Sequence_Names
#and subsequent lines like (tab-delimited):
#18S	Lewisia_cotyledon_cotyledon_14	1	Lewisia_cotyledon_cotyledon_14_1308_1249_14445_R1;Lewisia_cotyledon_cotyledon_14_1311_17462_24991_R1....
#[0]: locus name
#[1]: individual name
#[2]: number of BLAST hits these sequences had
#[3]: semi-colon-delimited list of sequence names (from fasta file)
#These sequence names will need to be converted to the actual sequence names in the BLAST file to find the sequences
#It writes one file with the sequences and one with the phred scores, as well
#as a script for the Minimo analysis.

import sys #to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels
import gzip #We want to be able to open zipped files.
from itertools import izip #We want to be able to look at two files simultaneously (without having to
#have them both fully in memory at the same time

#Example
#tblast_to_fasta_qual.py ~/transcriptomes/TS27_2015_01_26/TS27_MiSeq_Inds_Copy/o2_Seqs_to_Loci.txt ~/transcriptomes/TS27_2015_01_26/TS27_MiSeq_Inds/ t_ ~/transcriptomes/TS27_2015_01_26/TS27_MiSeq_Inds_Copy/blast_hits_Min/ b5_ separate qual
#tblast_to_fasta_qual.py InFileName SeqFolder SeqFilePre OutFolder OutFilePre OutMode [together, separate]

Usage ='''
tblast_to_fasta_qual.py version 1.0 makes new fasta files for individual loci
as well as separate files with the phred scores from the blast output parsed by 
tbaits_blastn_parse.py or tbaits_blast_parse.py, as well as a Minimo script
for analysis of these files:
tblast_to_fastq.py [infile name] [sequence folder] [prefix for sequence files]
[out folder] [out file prefix or "none" for none]
[mode: "together" if all reads for the same locus should be analyzed together, 
even if they come from different individuals or "separate" if reads for each locus
should be analyzed separately for the different individuals] '''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 8:
	sys.exit("ERROR!  This script requires 8 additional arguments and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
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
NCores = sys.argv[7]

######################################################################################

#OrderedListMaking makes a list in which objects are arranged in descending 
#(if Rev == True) or ascending (if Rev == False) order
#It takes a dictionary with DictTemp[Thing] = Size
#original
def OrderedListMaking(DictTemp, Rev):
	DictbySize = defaultdict(list)
	for Thing in DictTemp:
		Size = DictTemp[Thing]
		DictbySize[Size].append(Thing)
	ListTemp = [ ]#list of Things in descending order of size
	for Size in sorted(DictbySize.keys(), reverse = Rev):
		ListTemp += DictbySize[Size]
	return(ListTemp)

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

######################################################################################

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
			SeqName = "_".join(SeqName[-3:])#This is the version for normal files
			#SeqName = "_".join(SeqName[-4:-1])#This is the version for TS27_1
			SeqNameDict[Ind][SeqName] = Locus
InFile.close()

IndList = SeqNameDict.keys()
IndList = sorted(IndList)

print("Information on sequences from %d individuals was read from the file %s.\n" % (len(IndList), InFileName))
sys.stderr.write("Information on sequences from %d individuals was read from the file %s.\n" % (len(IndList), InFileName))

NumLocusSeqs = defaultdict(int)#NumLocusSeqs[Locus] = NumSeqs
#for each individual at once:
#read fastq files and make defaultdict indseqs[seqname]['R1' or 'R2'] = fastq sequence
#while doing this, check both sequences to see if they are entirely Ns, in which case
#write "error" for the fastq sequence
for Ind in IndList:
	IndSeqs = defaultdict(dict) #dictionary of sequences for that individual
	IndSeqs['R1'] = defaultdict(dict)
	IndSeqs['R2'] = defaultdict(dict)
	IndSeqs['S'] = defaultdict(dict)
	IndQuals = defaultdict(dict) #dictionary of quality scores for those sequences
	IndQuals['R1'] = defaultdict(dict)
	IndQuals['R2'] = defaultdict(dict)
	IndQuals['S'] = defaultdict(dict)
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
				NumLocusSeqs[Locus] += 1
			except KeyError:
				SeqWanted = False
		#check to see if the sequences are mainly Ns
		if (SeqLine == 1) and (SeqWanted == True):
			NumNs1 = Line1.count('N')
			if NumNs1 > 3:
				Seq1 = "error"
				Seq1Good = False
			else:
				#If not, add the sequence to the dictionary.
				Seq1 = Line1
				Seq1Good = True
			NumNs2 = Line2.count('N')
			if NumNs2 > 3:
				Seq2 = "error"
				Seq2Good = False
			else:
				Seq2 = Line2
				Seq2Good = True
		#skip the third line
		#add the fourth line
		elif (SeqLine == 3) and (SeqWanted == True):
			Line1P = list(Line1)
			Line1N = [ ]
			for Num in Line1P:
				Line1N.append(str(ord(Num)))
			Line1 = " ".join(Line1N)			
			Line2P = list(Line2)
			Line2N = [ ]
			for Num in Line2P:
				Line2N.append(str(ord(Num)))
			Line2 = " ".join(Line2N)
			#check to see if we want the sequence before adding it to the dictionary of sequences to be written
			if (Seq1Good == True) and (Seq2Good == True):
				IndSeqs['R1'][Locus][SeqName] = Seq1
				IndSeqs['R2'][Locus][SeqName] = Seq2
				IndQuals['R1'][Locus][SeqName] = Line1
				IndQuals['R2'][Locus][SeqName] = Line2
				BothGood += 1
			elif (Seq1Good == True) and (Seq2Good == False):
				IndSeqs['S'][Locus][SeqName] = Seq1
				IndQuals['S'][Locus][SeqName] = Line1
				OneGood += 1
			elif (Seq1Good == False) and (Seq2Good == True):
				IndSeqs['S'][Locus][SeqName] = Seq2
				IndQuals['S'][Locus][SeqName] = Line2
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
			OutFileName1 = OutFolder+OutFilePre+Locus+"_"+Ind+"_seqs.fa" 
			OutFile1 = open(OutFileName1, 'w')		
			OutFileName2 = OutFolder+OutFilePre+Locus+"_"+Ind+"_phred.txt" 
			OutFile2 = open(OutFileName2, 'w')
			OutScriptDict[Locus][Ind]=[OutFileName1, OutFileName2]
		elif OutMode == "together":
			OutFileName1 = OutFolder+OutFilePre+Locus+"_seqs.fa" 
			OutFile1 = open(OutFileName1, 'a')		
			OutFileName2 = OutFolder+OutFilePre+Locus+"_phred.txt" 
			OutFile2 = open(OutFileName2, 'a')
			OutScriptDict[Locus] = [OutFileName1, OutFileName2] 
		for SeqName in IndSeqs['R1'][Locus]:
			OutLine1 = ">"+Ind+"_"+SeqName+"\n"+IndSeqs['R1'][Locus][SeqName]+"\n"
			OutFile1.write(OutLine1)
			OutLine1 = ">"+Ind+"_"+SeqName+"\n"+IndSeqs['R2'][Locus][SeqName]+"\n"
			OutFile1.write(OutLine1)
			OutLine2 = ">"+Ind+"_"+SeqName+"\n"+IndQuals['R1'][Locus][SeqName]+"\n"
			OutFile2.write(OutLine2)
			OutLine2 = ">"+Ind+"_"+SeqName+"\n"+IndQuals['R2'][Locus][SeqName]+"\n"
			OutFile2.write(OutLine2)
		try:
			for SeqName in IndSeqs['S'][Locus]:
				OutLine1 = ">"+Ind+"_"+SeqName+"\n"+IndSeqs['S'][Locus][SeqName]+"\n"
				OutFile1.write(OutLine1)
				OutLine2 = ">"+Ind+"_"+SeqName+"\n"+IndQuals['S'][Locus][SeqName]+"\n"
				OutFile2.write(OutLine2)
		except KeyError:
			"do nothing"
		OutFile1.close()
		OutFile2.close()
		NumFiles += 2
	print("Sequences were written to %d fasta and quality score files, with \nnames such as %s and %s.\n" \
	% (NumFiles, OutFileName1, OutFileName2))
	sys.stderr.write("Sequences were written to %d fasta and quality score files, with \nnames such as %s and %s.\n" \
	% (NumFiles, OutFileName1, OutFileName2))

#Making a script to analyze the sequences using Minimo:
#first ordering the loci by the number of reads they have
OutLocusList = OrderedListMaking(NumLocusSeqs, True)
#then writing the script with the largest loci first
OutScript1 = [ "#! /bin/bash\n\n" ]
OutScript2 = [ ]
for Locus in OutLocusList:
	if OutMode == "separate":
		for Ind in OutScriptDict[Locus]:
			Line = "Minimo "+OutScriptDict[Locus][Ind][0]+" -D QUAL_IN="+OutScriptDict[Locus][Ind][1]+" -D FASTA_EXP=1 -D MIN_LEN=40 -D MIN_IDENT=90\n"
			OutScript2.append(Line)
	elif OutMode == "together":
		Line = "Minimo "+OutScriptDict[Locus][0]+" -D QUAL_IN="+OutScriptDict[Locus][1]+" -D FASTA_EXP=1 -D MIN_LEN=40 -D MIN_IDENT=90\n"
		OutScript2.append(Line)

OutFileName1 = OutFolder+"minimo_script_"+OutMode+"1.sh"
OutFileName2 = OutFolder+"minimo_script_"+OutMode+"2.sh"

Line = "cat "+OutFileName2+" | parallel --jobs "+NCores+" --joblog "+OutFolder+"parallel_log.log\n"
OutScript1.append(Line)

OutFileWriting(OutFileName1, OutScript1)
OutFileWriting(OutFileName2, OutScript2)

if OutMode == "separate":
	print("The shell script to analyze the individuals separately using Minimo is %s.\n" % (OutFileName1))
	sys.stderr.write("The shell script to analyze the individuals separately using Minimo is %s.\n" % (OutFileName1))
elif OutMode == "together":
	print("The shell script to analyze the all individuals together using Minimo is %s.\n" % (OutFileName1))
	sys.stderr.write("The shell script to analyze the all individuals together using Minimo is %s.\n" % (OutFileName1))

LocusList = OutScriptDict.keys()
LocusList = sorted(LocusList)

OutFileName = OutFolder+OutFilePre+"Locus_List.txt"
OutFile = open(OutFileName, 'w')
for Locus in LocusList:
	Line = Locus+"\t"+",".join(OutScriptDict[Locus].keys())+"\n"
	OutFile.write(Line)
OutFile.close()

print("The list of loci and the individuals that have sequences for those loci was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The list of loci and the individuals that have sequences for those loci was written to the file %s.\n" % (OutFileName))

