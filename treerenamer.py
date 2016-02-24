#! /usr/bin/env python

#treerenamer.py version 1.0 6 Oct. 2015 Abby Moore
#treerenamer.py takes a Newick-formatted tree file and changes the names.
#It simply looks for the old names and replaces them with the new names, so if there are other lines in the file
#that do not have any of the old names in them, it will just copy those unaltered.
#It just can't handle really large files, because it saves everything as a list before writing it to the output file.

import sys
import re

Usage='''
treerenamer.py version 1.0
This file will rename the taxa in a tree (or an alignment or any other file)
according to an input list of the format:
oldname[tab]newname
treerenamer.py
[input file with trees]
[name for output file]
[file with the old and new names]
'''

'''
treerenamer.py InFileName OutFileName NDictFileName
treerenamer.py ~/Documents/Molluginaceae/Mollugo201510/RAxML_bipartitions.PortITS ~/Documents/Molluginaceae/Mollugo201510/RAxML_bipartitions.PortITSnn ~/Documents/Molluginaceae/newnames_20151006.txt 
'''

if len(sys.argv) != 4:
	sys.exit("ERROR!!  This script takes 3 additional arguments and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFileName = sys.argv[1]
OutFileName = sys.argv[2]
NameDictFileName = sys.argv[3]

#DictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
#from tbaits_intron_removal.py
def DictFromFile(FileName):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]] = Line[1]
	InFile.close()
	#print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	#sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is NameDict

#OutFileWritingN writes an output file from a list of lines to write.
#This version adds "\n" to the ends of lines
#from OutFileWriting from tbaits_intron_removal.py
def OutFileWritingN(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line+"\n")
	OutFile.close()
	#print("Output file %s written.\n" % (FileName))
	sys.stderr.write("Output file %s written.\n" % (FileName))

#TreeRenamer renames a tree file (or an alignment or any file in which these names occur) according to
#the a dictionary of names NDict[OldName] = NewName
def TreeRenamer(TFName,NDict):
	#opening the tree file
	OList = [ ]
	InFile = open(TFName, 'rU')
	for Line in InFile:
		TreeString = Line.strip('\n').strip('\r')
		for Name in NDict:
			TreeString = TreeString.replace(Name,NDict[Name])
		OList.append(TreeString)
	InFile.close()
	return(OList)
	#This is OutList

NameDict = DictFromFile(NameDictFileName)
OutList =  TreeRenamer(InFileName, NameDict)
OutFileWritingN(OutFileName,OutList)
