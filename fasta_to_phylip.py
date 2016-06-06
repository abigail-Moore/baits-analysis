#! /usr/bin/env python

#fasta_to_phylip.py version 1.0 20 April 2015 Abby Moore
#This is just a very simple script to convert an alignment in fasta format to extended
#phylip format.  It assumes the name of the phylip file will be the same as that of
#the fasta file, but that the extensions will be different (.phy and .fa, respectively).

import sys
from Bio import AlignIO

Usage = '''
fasta_to_phylip.py is a script to convert a fasta alignment to a phylip
alignment (hopefully) without messing up the names too badly.
fasta_to_phylip.py [name of fasta alignment, expects an extension of .fa]
'''

if len(sys.argv) != 2:
	sys.exit("ERROR! This script expects one additional argument, and you gave it %d arguments!  %s" % (len(sys.argv), Usage))
InFileName = sys.argv[1]

MyAlignment = AlignIO.read(InFileName, "fasta")
for seq_record in MyAlignment:
	SeqNameTemp = seq_record.id
	for Char in SeqNameTemp:
		if (Char in [':',',','(',')',':','[',']',"'"]):
			SeqNameTemp = SeqNameTemp.replace(Char,'-')
	seq_record.id = SeqNameTemp

OutFileName = InFileName[:-2]+"phy"

AlignIO.write(MyAlignment,OutFileName,'phylip-relaxed')
