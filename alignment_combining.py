#! /usr/bin/env python

#alignment_combining.py version 1.0 13 Nov. 2016 Abby Moore
#This script combines two alignments in whatever format you desire, while removing any sequences with the same names.
#When duplicates are found, it keeps the sequence from the first alignment.

import sys
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

Usage = '''
alignment_combining.py version 1.0 13 Nov. 2016
This script combines two alignments (or sets of unaligned sequences), while 
removing sequences from the second alignment that are already present in the
first
alignment_combining.py
[name of first alignment]
[format of first alignment]
[name of second alignment]
[format of second alignment]
[name of output alignment]
[format of output alignment]
'''

#alignment_combining.py Al1Name Al1Format Al2Name Al2Format AlOutName AlOutFormat

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 7:
	sys.exit("ERROR!  seqsorting.py requires 6 additional arguments and you supplied %d.\n %s" % (len(sys.argv), Usage))
Al1Name = sys.argv[1]
Al1Format = sys.argv[2]
Al2Name = sys.argv[3]
Al2Format = sys.argv[4]
AlOutName = sys.argv[5]
AlOutFormat = sys.argv[6]

##########################################################################################

#SeqFiletoDict
#This reads a sequence file of a specified format and makes a dictionary of the sequences.
def SeqFiletoDict(FileName, FileFormat):
	DictTemp = { }
	InFile = open(FileName, 'rU')
	for record in SeqIO.parse(InFile, FileFormat):
		DictTemp[record.id] = str(record.seq)
	InFile.close()
	print("%d sequences were read from the file %s.\n" % (len(DictTemp.keys()), FileName))
	#sys.stderr.write("%d sequences were read from the file %s.\n" % (len(DictTemp.keys()), FileName)
	return DictTemp

#########################################################################################################

SeqDict1 = SeqFiletoDict(Al1Name, Al1Format)
SeqDict2 = SeqFiletoDict(Al2Name, Al2Format)

#add the sequences from SeqDict2 to SeqDict1, if a sequence with that name is not already present
for SeqName in SeqDict2:
	try:
		SeqTemp = SeqDict1[SeqName]
	except KeyError:
		SeqDict1[SeqName] = SeqDict2[SeqName]

#sorting the sequences:
SeqNameList = sorted(SeqDict1.keys())

#Writing the sequences:
SeqsWritten = 0
OutFile = open(AlOutName, 'w')
for SeqName in SeqNameList:
	Record1 = SeqRecord(seq=Seq(SeqDict1[SeqName], IUPAC), id = SeqName, description = "")
	SeqIO.write(Record1, OutFile, AlOutFormat)
	SeqsWritten += 1
OutFile.close()
print("%d sequences were written to the file %s.\n" % (SeqsWritten, AlOutName))
#sys.stderr.write("%d sequences were written to the file %s.\n" % (SeqsWritten, AlOutName))