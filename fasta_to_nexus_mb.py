#! /usr/bin/env python

#fasta_to_nexus_mb.py version 1.0 20 April 2015 Abby Moore
#This is just a very simple script to convert an alignment in fasta format to nexus
#format.  It assumes the name of the nexus file will be the same as that of
#the fasta file, but that the extensions will be different (.nex and .fa, respectively).
#It is the same as fasta_to_nexus.py, but it also makes a separate file (FileName_mb.nex)
#that will execute this file in Mr.Bayes.

import sys
from Bio import AlignIO
from Bio import SeqIO #to read files of non-aligned sequences
from Bio.Alphabet import generic_dna #to recognize sequences
from Bio.SeqRecord import SeqRecord #to make strings into sequence objects
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq #to edit sequences


Usage = '''
fasta_to_nexus_mb.py is a script to convert a fasta alignment to a nexus
alignment (hopefully) without messing up the names too badly.
It also makes a separate Mr.Bayes script..
fasta_to_nexus.py [name of fasta alignment, expects an extension of .fa]
'''

if len(sys.argv) != 2:
	sys.exit("ERROR! This script expects one additional argument, and you gave it %d arguments!  %s" % (len(sys.argv), Usage))
InFileName = sys.argv[1]

sys.stderr.write("Alignment %s will be processed.\n" % (InFileName))

MyAlignment = MultipleSeqAlignment([], generic_dna)
MyAlignment.extend(AlignIO.read(InFileName, 'fasta'))
for record in MyAlignment:
	SeqNameTemp = record.id
	for Char in SeqNameTemp:
		if (Char in [':',',','(',')',':','[',']',"'","="]):
			SeqNameTemp = SeqNameTemp.replace(Char,'-')
	record.id = SeqNameTemp

OutFileName = InFileName[:-2]+"nex"

AlignIO.write(MyAlignment,OutFileName,'nexus')

OutList = [ ]
Line = "set autoclose=yes nowarn=yes\nexecute "+OutFileName+"\n"
Line += "lset nst=6 rates=invgamma\nunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all)\n"
OutList.append(Line)
Line = "prset ratepr=variable\nmcmcp ngen= 4000000 relburnin=yes burninfrac=0.25  printfreq=10000  samplefreq=1000 nruns=2 nchains=4 savebrlens=yes\n"
OutList.append(Line)
Line = "mcmc\nsumt\nquit\n"
OutList.append(Line)

OutFileName = InFileName[:-3]+"_mb.nex"
OutFile = open(OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()
