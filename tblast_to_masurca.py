#! /usr/bin/env python

#tblast_to_masurca.py version 1.0 23 May 2015
#This script reads the Seqs_to_Loci.txt output 
#from tbaits_blastn_parse.py or tbaits_blast_parse.py
#and writes files for each locus for each individual.
#modified from tblast_to_fastq.py.
#It expects file with the header:
#Locus_Name	Individual_Name	Number_of_Hits	Sequence_Names
#and subsequent lines like (tab-delimited):
#18S	Lewisia_cotyledon_cotyledon_14	1	Lewisia_cotyledon_cotyledon_14_1308_1249_14445_R1;Lewisia_cotyledon_cotyledon_14_1311_17462_24991_R1....
#[0]: locus name
#[1]: individual name
#[2]: number of BLAST hits these sequences had
#[3]: semi-colon-delimited list of sequence names (from fasta file)
#These sequence names will need to be converted to the actual sequence names in the BLAST file to find the sequences
#It writes R1 and R2 fastq files and a separate config file for each Masurca analysis
#(one for each locus for each individual), as well as a bash script to run everything
#This only works in the "separate" mode

import sys #to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels
import gzip #We want to be able to open zipped files.
from itertools import izip #We want to be able to look at two files simultaneously (without having to
#have them both fully in memory at the same time)

#examples
'''
tblast_to_masurca.py InFileName SeqFolder SeqFilePre OutFolder OutFilePre
'''

Usage ='''
tblast_to_masurca.py version 1.0 makes new fastq files for individual loci from
the blast output parsed by tbaits_blastn_parse.py or tbaits_blast_parse.py and
masurca config files for each locus/individual
tblast_to_fastq.py [infile name] [sequence folder] [prefix for sequence files]
[out folder] [out file prefix or "none" for none]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 6:
	sys.exit("Error!!  tblast_to_fastq.py requires 5 additional arguments and you supplied %d.  %s" % (len(sys.argv)-1, Usage))
else:
	InFileName = sys.argv[1]
	SeqFolder = sys.argv[2]
	SeqFilePre = sys.argv[3]
	OutFolder = sys.argv[4]
	if sys.argv[5] == "none":
		OutFilePre = ""
	else:
		OutFilePre = sys.argv[5]

if SeqFolder[-1] != "/":
	SeqFolder += "/"

SeqNameDict = defaultdict(dict) #The dictionary that will have the sequence names from each locus
OutScriptDict = defaultdict(dict) #The dictionary that will be used to make the outfile.

#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
#from tbaits_intron_removal.py
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	#print("Output file %s written.\n" % (FileName))
	#sys.stderr.write("Output file %s written.\n" % (FileName))

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
			SeqName = "_".join(SeqName[-3:])
			SeqNameDict[Ind][SeqName] = Locus
InFile.close()

IndList = SeqNameDict.keys()
IndList = sorted(IndList)

print("Information on sequences from %d individuals was read from the file %s.\n" % (len(IndList), InFileName))
sys.stderr.write("Information on sequences from %d individuals was read from the file %s.\n" % (len(IndList), InFileName))

#for each individual at once:
#read fastq files and make defaultdict indseqs[seqname]['R1' or 'R2'] = fastq sequence
#while doing this, check both sequences to see if they are entirely Ns, in which case
#write "error" for the fastq sequence
for Ind in IndList:
	IndSeqs = defaultdict(dict) #dictionary of sequences for that individual
	IndSeqs['R1'] = defaultdict(dict)
	IndSeqs['R2'] = defaultdict(dict)
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
				Seq1 = Line1 + "\n"
				Seq2 = Line2 + "\n"
			except KeyError:
				SeqWanted = False
		#check to see if the sequences are mainly Ns
		if (SeqLine == 1) and (SeqWanted == True):
			NumNs1 = Line1.count('N')
			if NumNs1 > 3:
				Seq1 = "error"
				Seq1Good = False
			else:
				Seq1 += Line1 + "\n"
				Seq1Good = True
			NumNs2 = Line2.count('N')
			if NumNs2 > 3:
				Seq2 = "error"
				Seq2Good = False
			else:
				Seq2 += Line2 + "\n"
				Seq2Good = True
		#add the third line
		elif (SeqLine == 2) and (SeqWanted == True):
			Seq1 += Line1 + "\n"
			Seq2 += Line2 + "\n"
		#add the fourth line
		elif (SeqLine == 3) and (SeqWanted == True):
			Seq1 += Line1 + "\n"
			Seq2 += Line2 + "\n"
			#check to see if we want the sequence before adding it to the dictionary of sequences to be written
			if (Seq1Good == True) and (Seq2Good == True):
				IndSeqs['R1'][Locus][SeqName] = Seq1
				IndSeqs['R2'][Locus][SeqName] = Seq2
				BothGood += 1
			elif (Seq1Good == True) and (Seq2Good == False):
				OneGood += 1
			elif (Seq1Good == False) and (Seq2Good == True):
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
		OutFileName1 = OutFolder+OutFilePre+Locus+"_"+Ind+"_R1.fastq" 
		OutFile1 = open(OutFileName1, 'w')		
		OutFileName2 = OutFolder+OutFilePre+Locus+"_"+Ind+"_R2.fastq" 
		OutFile2 = open(OutFileName2, 'w')
		OutScriptDict[Locus][Ind]=[OutFilePre+Locus+"_"+Ind+"_R1.fastq", OutFilePre+Locus+"_"+Ind+"_R2.fastq"]
		for SeqName in IndSeqs['R1'][Locus]:
			OutFile1.write(IndSeqs['R1'][Locus][SeqName])
			OutFile2.write(IndSeqs['R2'][Locus][SeqName])
		OutFile1.close()
		OutFile2.close()
		NumFiles += 2
	print("Sequences were written to %d files, with names such as %s.\n" % (NumFiles, OutFileName1))
	sys.stderr.write("Sequences were written to %d files, with names such as %s.\n" % (NumFiles, OutFileName1))

ConfigScriptEnd = '''END

PARAMETERS
#this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for Illumina-only assemblies and to 0 if you have 1x or more long (Sanger, 454) reads, you can also set this to 0 for large data sets with high jumping clone coverage, e.g. >50x
USE_LINKING_MATES = 1
#this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS = cgwErrorRate=0.15 ovlMemory=4GB
#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if coverage >100
KMER_COUNT_THRESHOLD = 1
#auto-detected number of cpus to use
NUM_THREADS = 4
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage
JF_SIZE = 50000000
#this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
DO_HOMOPOLYMER_TRIM = 0
END
'''

#Making a script to analyze the sequences using Masurca and making the individual config files:
OutList = ["#! /bin/bash\n\n"]
Line = "mkdir "+OutFolder+OutFilePre+"Masurca_Contigs/\n"
OutList.append(Line)
for Locus in OutScriptDict.keys():
	Line = "mkdir "+OutFolder+Locus+"\n"
	OutList.append(Line)
	for Ind in OutScriptDict[Locus]:
		OF2 = OutFolder+Locus+"/"+Locus+"_"+Ind+"/"
		Line = "mkdir "+OF2+"\n"
		Line += "mv "+OutFolder+OutScriptDict[Locus][Ind][0]+" "+OF2+"\n"
		Line += "mv "+OutFolder+OutScriptDict[Locus][Ind][1]+" "+OF2+"\n"
		Line += "mv "+OutFolder+OutFilePre+Locus+"_"+Ind+"_config.sh "+OF2+"\n"
		Line += "cd "+OF2+"\n"
		Line += "masurca "+OutFilePre+Locus+"_"+Ind+"_config.sh\n"
		Line += "./assemble.sh\n"
		OutList.append(Line)
		ConfigScript = [ 'DATA\n']
		ConfigScript.append('PE= pe 180 25  '+OF2+OutScriptDict[Locus][Ind][0]+" "+OF2+OutScriptDict[Locus][Ind][1]+" "+"\n")
		ConfigScript.append(ConfigScriptEnd)
		OutFileName = OutFolder+OutFilePre+Locus+"_"+Ind+"_config.sh"
		OutFileWriting(OutFileName, ConfigScript)

#writing the bash script:
OutFileName = OutFolder+"masurca_script.sh"
OutFileWriting(OutFileName, OutList)
print("The shell script to analyze the sequences with masurca is %s.\n" % (OutFileName))
sys.stderr.write("The shell script to analyze the sequences with masurca is %s.\n" % (OutFileName))

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

''' Goal:
# example configuration file 

# DATA is specified as type {PE,JUMP,OTHER} and 5 fields:
# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
# innies, i.e. --->.<---, and JUMP are assumed to be outties
# <---.--->. If there are any jump libraries that are innies, such as
# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
# are optional for PE libraries and mandatory for JUMP libraries. Any
# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
# converted into Celera Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
PE= pe 180 25  /users/ajm3/scratch/eedwards/TS31_1/blast_gfams/masurca_tr_ppc/bf1_ppc_Lewisia_cotyledon_cotyledon_14_R1.fastq /users/ajm3/scratch/eedwards/TS31_1/blast_gfams/masurca_tr_ppc/bf1_ppc_Lewisia_cotyledon_cotyledon_14_R2.fastq
END

PARAMETERS
#this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for Illumina-only assemblies and to 0 if you have 1x or more long (Sanger, 454) reads, you can also set this to 0 for large data sets with high jumping clone coverage, e.g. >50x
USE_LINKING_MATES = 0
#this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS = cgwErrorRate=0.15 ovlMemory=4GB
#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if coverage >100
KMER_COUNT_THRESHOLD = 1
#auto-detected number of cpus to use
NUM_THREADS = 4
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage
JF_SIZE = 500000
#this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
DO_HOMOPOLYMER_TRIM = 0
END
'''
