#! /usr/bin/env python

#seqsorting.py version 1.0 3 July 2015 Abby Moore
#This script takes an alignment in whatever format you desire and sorts the sequences alphabetically.

import sys
from Bio.Seq import Seq #to edit sequences
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import IUPAC #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects

Usage = '''
seqsorting.py version 1.0
This script sorts the sequences in a sequence file alphabetically.
seqsorting.py
[sequence format--in]
[sequence format--out, or same if the same is the input format]
[infile name]
[outfile name]
'''

#seqsorting.py InFormat OutFormat InFileName OutFileName 

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 5:
	sys.exit("ERROR!  seqsorting.py requires 4 additional arguments and you supplied %d.\n %s" % (len(sys.argv), Usage))
InFormat = sys.argv[1]
OutFormat = sys.argv[2]
InFileName = sys.argv[3]
OutFileName = sys.argv[4]

SeqDict = { }

#Reading the sequences:
InFile = open(InFileName, 'rU')
for record in SeqIO.parse(InFile, InFormat):
	SeqDict[record.id] = str(record.seq)
InFile.close()
print("%d sequences were read from the file %s.\n" % (len(SeqDict.keys()), InFileName))
#sys.stderr.write("%d sequences were read from the file %s.\n" % (len(SeqDict.keys()), InFileName))

#sorting the sequences:
SeqNameList = sorted(SeqDict.keys())

#Writing the sequences:
SeqsWritten = 0
OutFile = open(OutFileName, 'w')
for SeqName in SeqNameList:
	Record1 = SeqRecord(seq=Seq(SeqDict[SeqName], IUPAC), id = SeqName, description = "")
	SeqIO.write(Record1, OutFile, OutFormat)
	SeqsWritten += 1
OutFile.close()
print("%d sequences were written to the file %s.\n" % (SeqsWritten, OutFileName))
#sys.stderr.write("%d sequences were written to the file %s.\n" % (SeqsWritten, OutFileName))
