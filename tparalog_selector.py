#! /usr/bin/env python

#tparalog_selector.py version 1.0 31 March 2016 Abby Moore
#This script goes through the output from tnotung_homolog_parsing.py and makes various lists of the loci
#and individuals to be used, depending on what percentage of the genes need to be sampled.

#File Formats:
#InFilePre_seqs_per_ind.txt
'''
Header Row: Individual_Name\tTotal_Sequences
IndName [0]: Alluaudia_dumosa_65
NumSeqs [1]: 218
'''
#InFilePre_Inds_per_Paralog.txt
'''
Header Row: Group\tloci
Header Column: Group\ngroup names\nTotal_Sequences_in_All_Groups\nNumber_of_Groups_with_Sequences
Group/category name [0]: Anacampserotaceae or Number_of_Groups_with_Sequences
NumSeqs [1]: 0 or 1
'''

import sys
from collections import defaultdict
import numpy #calculating means and percentiles

'''
tparalog_selector.py InFolder InFilePre OutFolder OutFilePre
'''

Usage = '''
tparalog_selector.py
This script selects the individuals and paralogs that will be used to construct
the final trees, according to various criteria.  It uses the output from
tnotung_homolog_parsing.py
[folder in which the output for tnotung_homolog_parsing.py is found]
[prefix for these files, or "none", if none]
[folder to which the files should be written]
[prefix for the output files]
'''

if len(sys.argv) < 5:
	sys.exit("ERROR!!!  This script requires at least 4 additional arguments, and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
if InFolder[-1] != "/":
	InFolder += "/"
InFilePre = sys.argv[2]
if InFilePre == "none":
	InFilePre = ""
OutFolder = sys.argv[3]
if OutFolder == "same":
	OutFolder = InFolder
elif OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[4]
if OutFilePre == "none":
	OutFilePre = ""

print(" ".join(sys.argv))

Verbose = False

#Verbose = True

#################################################################################

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

#################################################################################

#reading the file showing how many sequences each individual had
InFileName = InFolder+InFilePre+"seqs_per_ind.txt"
NumIndSeqs = { }
NumIndSeqList = [ ]
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	if Line[0] !=  'Individual_Name':
		NumIndSeqs[Line[0]] = int(Line[1])
		NumIndSeqList.append(int(Line[1]))
InFile.close()
print("File %s was read.\n" % (InFileName))
sys.stderr.write("File %s was read.\n" % (InFileName))

#then figuring out the mean and percentiles for that file
MeanNumSeqs = numpy.mean(NumIndSeqList)
Perc75 = numpy.percentile(NumIndSeqList, 75)
Perc25 = numpy.percentile(NumIndSeqList, 25)
Perc10 = numpy.percentile(NumIndSeqList, 10)
print("The mean number of sequences per individual is %.1f, and the interquartile range goes from %.1f to %.1f.\n" % (MeanNumSeqs, Perc25, Perc75))
sys.stderr.write("The mean number of sequences per individual is %.1f, and the interquartile range goes from %.1f to %.1f.\n" % (MeanNumSeqs, Perc25, Perc75))
MeanMinus2IQR = MeanNumSeqs-2*(Perc75-Perc25)
MeanPlus2IQR = MeanNumSeqs+2*(Perc75-Perc25)
if Verbose == True: print("Outliers will have more than %.1f or fewer than %.1f sequences.\n" % (MeanPlus2IQR, MeanMinus2IQR))
#figuring out which sequences are outliers or have below a certain percent of the total sequences
UpperOLList = [ ]
LowerOLList = [ ]
AbovePerc25 = [ ]
BelowPerc25 = [ ]
AbovePerc10 = [ ]
BelowPerc10 = [ ]
NonLowerOutlier = [ ]
for Ind in NumIndSeqs:
	#if an individual is a lower outlier, add it to the list
	if NumIndSeqs[Ind] < MeanMinus2IQR:
		LowerOLList.append(Ind)
	#if not, add it to the other list
	else:
		NonLowerOutlier.append(Ind)
		#and check to see if it is an upper outlier
		if NumIndSeqs[Ind] > MeanPlus2IQR:
			UpperOLList.append(Ind)
		#and check to see if it is above the 25th percentile
	if NumIndSeqs[Ind] > Perc25:
		AbovePerc25.append(Ind)
	else:
		BelowPerc25.append(Ind)
	#also check to see if it is above the 10th percentile
	if NumIndSeqs[Ind] > Perc10:
		AbovePerc10.append(Ind)
	else:
		BelowPerc10.append(Ind)
if Verbose == True: print("The following %d individuals were outliers on the high end: %s.\n" % (len(UpperOLList), ", ".join(UpperOLList)))
if Verbose == True: print("The following %d individuals were outliers on the low end: %s.\n" % (len(LowerOLList), ", ".join(LowerOLList)))
print("Of the %d individuals, %d fell above the 25th perctile, while %d fell above the 10th percentile, and %d were not lower outliers.\n" % (len(NumIndSeqs), len(AbovePerc25), len(AbovePerc10), len(NonLowerOutlier))) 

if len(NumIndSeqs) != (len(AbovePerc25)+len(BelowPerc25)):
	sys.exit("ERRORR!! There are missing individuals that have been classified as neither above nor below the 25th percentile!")
if len(NumIndSeqs) != (len(AbovePerc10)+len(BelowPerc10)):
	sys.exit("ERRORR!! There are missing individuals that have been classified as neither above nor below the 10th percentile!")

#printing the lists of individuals to files
OutFileName = OutFolder+OutFilePre+"above_25.txt"
OutList = "\n".join(AbovePerc25)
OutFileWriting(OutFileName, OutList)
OutFileName = OutFolder+OutFilePre+"above_10.txt"
OutList = "\n".join(AbovePerc10)
OutFileWriting(OutFileName, OutList)
OutFileName = OutFolder+OutFilePre+"below_25.txt"
OutList = "\n".join(BelowPerc25)
OutFileWriting(OutFileName, OutList)
OutFileName = OutFolder+OutFilePre+"below_10.txt"
OutList = "\n".join(BelowPerc10)
OutFileWriting(OutFileName, OutList)
OutFileName = OutFolder+OutFilePre+"non_lower_outliers.txt"
OutList = "\n".join(NonLowerOutlier)
OutFileWriting(OutFileName, OutList)
OutFileName = OutFolder+OutFilePre+"lower_outliers.txt"
OutList = "\n".join(LowerOLList)
OutFileWriting(OutFileName, OutList)

print("Six output files were written, with names such as %s, with the lists of individuals that conform (or do not conform) to various criteria.\n" % (OutFileName))
sys.stderr.write("Six output files were written, with names such as %s, with the lists of individuals that conform (or do not conform) to various criteria.\n" % (OutFileName))


#now looking at how many individuals or groups have each locus
InFileName = InFolder+InFilePre+"Inds_per_Paralog.txt"
#ParGroupDict/ParGroupList are the number of groups each paralog has sequences for
#ParIndDict/ParIndList are the number of individuals each paralog has sequences for
ParGroupDict = { }
ParGroupList = [ ]
ParIndDict = { }
ParIndList = [ ]
InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	if Line[0] == "Group":
		KeyLine = Line
	elif Line[0] == "Total_Sequences_in_All_Groups":
		for LinePos in range(1,len(Line)):
			NumInds = int(Line[LinePos])
			ParIndList.append(NumInds)
			ParIndDict[KeyLine[LinePos]] = NumInds
	elif Line[0] == "Number_of_Groups_with_Sequences":
		for LinePos in range(1,len(Line)):
			NumGroups = int(Line[LinePos])
			ParGroupList.append(NumGroups)
			ParGroupDict[KeyLine[LinePos]] = NumGroups
InFile.close()

PIMean = numpy.mean(ParIndList)
PIPerc25 = numpy.percentile(ParIndList, 25)
PIPerc10 = numpy.percentile(ParIndList, 10)
print("The mean number of individuals in which a paralog is found is %.2f, 75 percent of paralogs had at least %d individuals, and 90 percent of paralogs had at least %d individuals.\n" % (PIMean, PIPerc25, PIPerc10))
if PIPerc25 == PIPerc10:
	print("The groups based on the number of individuals with a given paralog will be the same when the 25th percentile or the 10th percentile is used as the cutoff value.\n")

PGMean = numpy.mean(ParGroupList)
PGPerc25 = numpy.percentile(ParGroupList, 25)
PGPerc10 = numpy.percentile(ParGroupList, 10)
print("The mean number of groups in which a paralog is found is %.2f, 75 percent of paralogs had at least %d groups, and 90 percent of paralogs had at least %d groups.\n" % (PGMean, PGPerc25, PGPerc10))
if PGPerc25 == PGPerc10:
	print("The groups based on the number of individuals with a given paralog will be the same when the 25th percentile or the 10th percentile is used as the cutoff value.\n")

#percentiles based on individuals with a given paralog
PIAbovePerc25 = [ ]
PIBelowPerc25 = [ ]
PIAbovePerc10 = [ ]
PIBelowPerc10 = [ ]
for ParName in ParIndDict:
	#check to see if it is above or below the 25th percentile
	if ParIndDict[ParName] > PIPerc25:
		PIAbovePerc25.append(ParName)
	else:
		PIBelowPerc25.append(ParName)
	#check to see if it is above or below the 10th percentile
	if ParIndDict[ParName] > PIPerc10:
		PIAbovePerc10.append(ParName)
	else:
		PIBelowPerc10.append(ParName)
#checking
if len(ParIndList) != (len(PIAbovePerc25)+len(PIBelowPerc25)):
	sys.exit("ERRORR!! There are missing loci that have been classified as neither above nor below the 25th percentile of number of individuals!")
if len(ParIndList) != (len(PIAbovePerc10)+len(PIBelowPerc10)):
	sys.exit("ERRORR!! There are missing loci that have been classified as neither above nor below the 10th percentile of number of individuals!")
#printing the lists
if Verbose == True: print("%d loci were present in more than %d individuals, the cutoff for the 25th percentile.\n" % (len(PIAbovePerc25), PIPerc25))
OutList = "\n".join(PIAbovePerc25)
OutFileName = OutFolder+OutFilePre+"Loci_above_25_ind.txt"
OutFileWriting(OutFileName, OutList)
if Verbose == True: print("%d loci were present in %d individuals or fewer, the cutoff for the 25th percentile.\n" % (len(PIBelowPerc25), PIPerc25))
OutList = "\n".join(PIBelowPerc25)
OutFileName = OutFolder+OutFilePre+"Loci_below_25_ind.txt"
OutFileWriting(OutFileName, OutList)
if Verbose == True: print("%d loci were present in more than %d individuals, the cutoff for the 10th percentile.\n" % (len(PIAbovePerc10), PIPerc10))
OutList = "\n".join(PIAbovePerc10)
OutFileName = OutFolder+OutFilePre+"Loci_above_10_ind.txt"
OutFileWriting(OutFileName, OutList)
if Verbose == True: print("%d loci were present in %d individuals or fewer, the cutoff for the 10th percentile.\n" % (len(PIBelowPerc10), PIPerc10))
OutList = "\n".join(PIBelowPerc10)
OutFileName = OutFolder+OutFilePre+"Loci_below_10_ind.txt"
OutFileWriting(OutFileName, OutList)

#percentiles based on groups with a given paralog
PGAbovePerc25 = [ ]
PGBelowPerc25 = [ ]
PGAbovePerc10 = [ ]
PGBelowPerc10 = [ ]
for ParName in ParGroupDict:
	#check to see if it is above or below the 25th percentile
	if ParGroupDict[ParName] > PGPerc25:
		PGAbovePerc25.append(ParName)
	else:
		PGBelowPerc25.append(ParName)
	#check to see if it is above or below the 10th percentile
	if ParGroupDict[ParName] > PGPerc10:
		PGAbovePerc10.append(ParName)
	else:
		PGBelowPerc10.append(ParName)
#checking
if len(ParIndList) != (len(PGAbovePerc25)+len(PGBelowPerc25)):
	sys.exit("ERRORR!! There are missing loci that have been classified as neither above nor below the 25th percentile of number of groups!")
if len(ParIndList) != (len(PGAbovePerc10)+len(PGBelowPerc10)):
	sys.exit("ERRORR!! There are missing loci that have been classified as neither above nor below the 10th percentile of number of groups!")
#printing the lists
if Verbose == True: print("%d loci were present in more than %d groups, the cutoff for the 25th percentile.\n" % (len(PGAbovePerc25), PGPerc25))
OutList = "\n".join(PGAbovePerc25)
OutFileName = OutFolder+OutFilePre+"Loci_above_25_group.txt"
OutFileWriting(OutFileName, OutList)
if Verbose == True: print("%d loci were present in %d groups or fewer, the cutoff for the 25th percentile.\n" % (len(PGBelowPerc25), PGPerc25))
OutList = "\n".join(PGBelowPerc25)
OutFileName = OutFolder+OutFilePre+"Loci_below_25_group.txt"
OutFileWriting(OutFileName, OutList)
if Verbose == True: print("%d loci were present in more than %d groups, the cutoff for the 10th percentile.\n" % (len(PGAbovePerc10), PGPerc10))
OutList = "\n".join(PGAbovePerc10)
OutFileName = OutFolder+OutFilePre+"Loci_above_10_group.txt"
OutFileWriting(OutFileName, OutList)
if Verbose == True: print("%d loci were present in %d groups or fewer, the cutoff for the 10th percentile.\n" % (len(PGBelowPerc10), PGPerc10))
OutList = "\n".join(PGBelowPerc10)
OutFileName = OutFolder+OutFilePre+"Loci_below_10_group.txt"
OutFileWriting(OutFileName, OutList)
print("Six output files were written, with names such as %s, with the lists of loci that conform (or do not conform) to various criteria.\n" % (OutFileName))
sys.stderr.write("Six output files were written, with names such as %s, with the lists of loci that conform (or do not conform) to various criteria.\n" % (OutFileName))
