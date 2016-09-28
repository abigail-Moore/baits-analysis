#! /usr/bin/env python

#tgenefam_to_spptree_master.py version 1.0 22 Feb. 2016 Abby Moore
#This is third master script after tfastq_assembly_master.py and tcontig_classif_master.py
#It goes from separate alignments for the various gene families that include both
#the sequences we want to analyze and various backbone sequences to species trees.

#depends on the following scripts:
'''
tcombpars_to_trees.py
tnotung_homolog_parsing.py
tfasta_to_phylip.py

and indirectly on:
mafft
raxml
notung
'''

import sys


Usage = '''
tgenefam_to_spptree_master.py runs all of the scripts associated with making 
species trees from the gene family alignments.
[the date to appear in the log file name]
[the folder in which the scripts are found, or none, if none]
[the path to the input/output folders]
[the main analysis prefix for the input/output folders]
[the suffix for the input alignments (from tparcomb_final.py)]
[the suffix for the output alignments]
[the file name for the list of loci]
[the file name for the file containing GroupName [tab] IndName]
[number of cores for Oscar]
[the file name of the species tree]
[the outgroup name, for use in making a rooted species tree-**note that there
must be an outgroup name for making a rooted spp tree if you are running in
mode New, or Notung will not work in the second round]
[mode, either Parallel or Array]
[tree mode, choose from New: make a new species tree using raxml analysis of the
concatenated sequences; Pause: pause after making gene trees, so that another
type of species tree can be made; Good: a good species tree was input, and 
notung does not need to be run for two rounds]
'''

'''
tgenefam_to_spptree_master.py Date ScriptFolder InOutFolder InOutFilePre AlFileInPost AlFileOutPost LLFileName IndGroupFileName NCores SpTrFileName OGName Mode[Parallel, Array] TreeMode[New, Pause, Good]
'''

ModeList = ['Parallel', 'Array']
TreeModeList = ['New', 'Pause', 'Good']

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 14:
	sys.exit("ERROR!  This script requires 13 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
Date = sys.argv[1]
ScriptFolder = sys.argv[2]
if ScriptFolder[-1] != "/":
	ScriptFolder += "/"
if ScriptFolder == "none/":
	ScriptFolder = ""
	ScriptFolder2 = "none"
else:
	ScriptFolder2 = ScriptFolder
InOutFolder = sys.argv[3]
if InOutFolder[-1] != "/":
	InOutFolder += "/"
InOutFilePre = sys.argv[4]
AlFileInPost = sys.argv[5]
AlFileOutPost = sys.argv[6]
LLFileName = sys.argv[7]
IndGroupFileName = sys.argv[8]
NCores = sys.argv[9]
SpTrFileName = sys.argv[10]
OGName = sys.argv[11]
Mode = sys.argv[12]
if Mode not in ModeList:
	sys.exit("ERROR!!  You wanted the mode %s, but the mode must be one of the following: %s\n%s" % (Mode, ", ".join(ModeList)))
TreeMode = sys.argv[13]
if TreeMode not in TreeModeList:
	sys.exit("ERROR!!  You wantede tree mode %s, but the tree mode must be one of the following: %s\n%s" % (TreeMode, ", ".join(TreeModeList)))

################################################################################

#NoneDictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value, but if the value is "none",
#then the value is ""
#modified from DictFromFile from tbaits_intron_removal.py
#not currently being used
def NoneDictFromFile(FileName):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]] = Line[1]
		if Line[1] == "none":
			TempDict[Line[0]] = ""
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is GroupDict

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
#not currently being used
def CaptureColumn(FileName, ColNum):
	TempList = [ ]
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r').split('\t')
		TempList.append(Line[ColNum])
	InFile.close()
	print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is PolyList

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

###################################################################################

if Mode == 'Parallel':
	OutFileName = InOutFolder+InOutFilePre+"_Species_Tree_Making.sh"
	OutScript = ["#! /bin/bash\n#SBATCH -J "+InOutFilePre+"_genefams\n#SBATCH -t 60:00:00\n#SBATCH -n "+NCores+"\n#SBATCH --mem="+str(int(NCores)*8)+"G\n"]
elif Mode == 'Array':
	OutScript = ["#! /bin/bash\n#SBATCH -J "+InOutFilePre+"_genefams1\n#SBATCH -t 2:00:00\n#SBATCH -n 1\n\n"]
	OutFileName = InOutFolder+InOutFilePre+"_Species_Tree_Making1.sh"
Line = "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\nmodule load mafft\nmodule load raxml\nmodule load mrbayes\n"
OutScript.append(Line)

##########tcombpars_to_trees.py
Line = "mkdir "+InOutFolder+InOutFilePre+"_genefams1\n"
if Mode == 'Parallel':
	Line += ScriptFolder+"tcombpars_to_trees.py "+LLFileName+" "+SpTrFileName+" "+InOutFolder+InOutFilePre+"_final/ "+InOutFilePre+"fi_ "+AlFileInPost+" "+InOutFolder+InOutFilePre+"_final2/ "+InOutFilePre+"fi2_ "+AlFileInPost+" "+InOutFolder+InOutFilePre+"_genefams1/ "+InOutFilePre+"gf1_ "+AlFileOutPost+" "+ScriptFolder2+" "+Mode+" >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	Line += "chmod u+x "+InOutFolder+InOutFilePre+"_genefams1/"+InOutFilePre+"gf1_analysis_script1.sh\n"
	Line += InOutFolder+InOutFilePre+"_genefams1/"+InOutFilePre+"gf1_analysis_script1.sh\n"
	Line += "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
elif Mode == 'Array':
	Line += ScriptFolder+"tcombpars_to_trees.py "+LLFileName+" "+SpTrFileName+" "+InOutFolder+InOutFilePre+"_final/ "+InOutFilePre+"fi_ "+AlFileInPost+" "+InOutFolder+InOutFilePre+"_final2/ "+InOutFilePre+"fi2_ "+AlFileInPost+" "+InOutFolder+InOutFilePre+"_genefams1/ "+InOutFilePre+"gf1_ "+AlFileOutPost+" "+ScriptFolder2+" "+Mode+" "+InOutFolder+InOutFilePre+"_Species_Tree_Making2.sh >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
	OutFileWriting(OutFileName, OutScript)
	OutFileName = InOutFolder+InOutFilePre+"_Species_Tree_Making2.sh"
	OutScript = ["#! /bin/bash\n#SBATCH -J "+InOutFilePre+"_genefams2\n#SBATCH -t 32:00:00\n#SBATCH -n "+NCores+"\n#SBATCH --mem="+str(int(NCores)*8)+"G\n\n"]
	Line = "module load mafft\nmodule load raxml\n"
	Line += "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
	
###########tnotung_homolog_parsing.py
Line = "mkdir "+InOutFolder+InOutFilePre+"_genetrees1\n"
Line += "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
Line += ScriptFolder+"tnotung_homolog_parsing.py "+InOutFolder+InOutFilePre+"_genefams1/ "+InOutFilePre+"gf1_ "+AlFileOutPost+" same "+InOutFilePre+"gf1_ "+AlFileOutPost+" "+LLFileName+" "+IndGroupFileName+" "+InOutFolder+InOutFilePre+"_genetrees1/ "+InOutFilePre+"gt1_ "+AlFileOutPost+" >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
Line += "chmod u+x "+InOutFolder+InOutFilePre+"_genetrees1/"+InOutFilePre+"gt1_analysis_script1.sh\n"
Line += InOutFolder+InOutFilePre+"_genetrees1/"+InOutFilePre+"gt1_analysis_script1.sh\n"
OutScript.append(Line)

###########tparalog_selector.py
Line = "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
Line += ScriptFolder+"tparalog_selector.py "+InOutFolder+InOutFilePre+"_genetrees1/ "+InOutFilePre+"gt1_ same "+InOutFilePre+"gt1_ >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
OutScript.append(Line)


#at this point, one of three things can happen:
Line = "mkdir "+InOutFolder+InOutFilePre+"_prunedgenetrees1\n"
OutScript.append(Line)
#1) if we don't have a good species tree and we want to run this automatically, we can make the concatenated tree and go on with that.
if TreeMode == "New":
	###########tal_combiner.py
	Line = ScriptFolder+"tal_combiner.py "+InOutFolder+InOutFilePre+"_genetrees1/"+InOutFilePre+"gt1_Loci_above_25_group.txt "+ScriptFolder+" "+InOutFolder+InOutFilePre+"_genetrees1/ "+InOutFilePre+"gt1_ "+AlFileOutPost+"_al.fa fasta fasta "+OGName+" "+InOutFolder+InOutFilePre+"_prunedgenetrees1 "+InOutFilePre+"pgt1_ combinednobs 10000 10000 all 2 Parallel >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	Line += "chmod u+x "+InOutFolder+InOutFilePre+"_prunedgenetrees1/"+InOutFilePre+"pgt1_analysis_script.sh\n"
	Line += InOutFolder+InOutFilePre+"_prunedgenetrees1/"+InOutFilePre+"pgt1_analysis_script.sh\n"
	Line += "cp "+InOutFolder+InOutFilePre+"_prunedgenetrees1/RAxML_bipartitions."+InOutFilePre+"pgt1_combined* "+InOutFolder+InOutFilePre+"pgt1_spptree.tre\n"
	OutScript.append(Line)
	#Then we need to go on to rerun tnotung_homolog_parsing.py, tparalog_selector.py, and tal_combiner.py again, but that overlaps with option #2, so this will be below.
#2) if we don't have a good species tree and we want to make one with astral or some other script oscar doesn't have, we can stop the script after the gene trees for astral have been made
elif TreeMode == "Pause":
	###########tal_combiner.py
	Line = ScriptFolder+"tal_combiner.py "+InOutFolder+InOutFilePre+"_genetrees1/"+InOutFilePre+"gt1_Loci_above_25_group.txt "+ScriptFolder+" "+InOutFolder+InOutFilePre+"_genetrees1/ "+InOutFilePre+"gt1_ "+AlFileOutPost+"_al.fa fasta fasta tree "+InOutFolder+InOutFilePre+"_prunedgenetrees1 "+InOutFilePre+"pgt1_ separate 10000 10000 all 2 "+Mode+" "+SpTrFileName+" >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
	if Mode == "Parallel":
		Line += "chmod u+x "+InOutFolder+InOutFilePre+"_prunedgenetrees1/"+InOutFilePre+"pgt1_analysis_script.sh\n"
		Line += InOutFolder+InOutFilePre+"_prunedgenetrees1/"+InOutFilePre+"pgt1_analysis_script.sh\n"
		OutScript.append(Line)
	OutFileWriting(OutFileName, OutScript)
	OutFileName = InOutFolder+InOutFilePre+"_New_Backbone_Spp_Tree.sh"
	print("After the Newick-formatted species tree has been made, save it as %s and run the script %s.\n" % (InOutFolder+InOutFilePre+"pgt1_spptree.tre", OutFileName))
	OutScript = ["#! /bin/bash\n#SBATCH -J "+InOutFilePre+"_genefams3\n#SBATCH -t 32:00:00\n#SBATCH -n "+NCores+"\n#SBATCH --mem="+str(int(NCores)*8)+"G\n\n"]
	Line = "module load mafft\nmodule load raxml\n"
	Line += "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
#3) if we do have a good species tree, then we don't need to restart the script after making one.
elif TreeMode == "Good":
	###########tal_combiner.py
	Line = ScriptFolder+"tal_combiner.py "+InOutFolder+InOutFilePre+"_genetrees1/"+InOutFilePre+"gt1_Loci_above_25_group.txt "+ScriptFolder+" "+InOutFolder+InOutFilePre+"_genetrees1/ "+InOutFilePre+"gt1_ "+AlFileOutPost+"_al.fa fasta fasta tree "+InOutFolder+InOutFilePre+"_prunedgenetrees1 "+InOutFilePre+"pgt1_ separate 10000 10000 "+InOutFolder+InOutFilePre+"_genetrees1/"+InOutFilePre+"gt1_non_lower_outliers.txt 2 "+Mode+" "+SpTrFileName+" >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
	Line += "chmod u+x "+InOutFolder+InOutFilePre+"_prunedgenetrees1/"+InOutFilePre+"pgt1_analysis_script.sh\n"
	Line += InOutFolder+InOutFilePre+"_prunedgenetrees1/"+InOutFilePre+"pgt1_analysis_script.sh\n"
	OutScript.append(Line)
	OutFileWriting(OutFileName, OutScript)

###########make concatenated tree or pause for species tree to be input
if (TreeMode == "New") or (TreeMode == "Pause"):
	#first I need to rerun Notung.  :(
	Line = "mkdir "+InOutFolder+InOutFilePre+"_genefams2\n"
	Line += "cp "+InOutFolder+InOutFilePre+"_genefams1/RAxML_bipartitions."+InOutFilePre+"gf1_*"+AlFileOutPost+" "+InOutFolder+InOutFilePre+"_genefams2/\n"
	Line += "ls "+InOutFolder+InOutFilePre+"_genefams2/RAxML_bipartitions."+InOutFilePre+"gf1_*"+AlFileOutPost+" | parallel java -jar "+ScriptFolder+"Notung-2.8.1.6-beta.jar -g {} -s "+InOutFolder+InOutFilePre+"pgt1_spptree.tre --speciestag prefix --rearrange --threshold 90 --nolosses --silent  --usegenedir\n"
	OutScript.append(Line)
	###########tnotung_homolog_parsing.py
	Line = "mkdir "+InOutFolder+InOutFilePre+"_genetrees2\n"
	Line += "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	Line += ScriptFolder+"tnotung_homolog_parsing.py "+InOutFolder+InOutFilePre+"_genefams2/ "+InOutFilePre+"gf1_ "+AlFileOutPost+" "+InOutFolder+InOutFilePre+"_genefams1/ "+InOutFilePre+"gf1_ "+AlFileOutPost+" "+LLFileName+" "+IndGroupFileName+" "+InOutFolder+InOutFilePre+"_genetrees2/ "+InOutFilePre+"gt2_ "+AlFileOutPost+" >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	Line += "chmod u+x "+InOutFolder+InOutFilePre+"_genetrees2/"+InOutFilePre+"gt2_analysis_script1.sh\n"
	Line += InOutFolder+InOutFilePre+"_genetrees2/"+InOutFilePre+"gt2_analysis_script1.sh\n"
	OutScript.append(Line)	
	###########tparalog_selector.py
	Line = "date >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	Line += ScriptFolder+"tparalog_selector.py "+InOutFolder+InOutFilePre+"_genetrees2/ "+InOutFilePre+"gt2_ same "+InOutFilePre+"gt2_ >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
	###########tal_combiner.py
	Line = "mkdir "+InOutFolder+InOutFilePre+"_prunedgenetrees2\n"
	Line += ScriptFolder+"tal_combiner.py "+InOutFolder+InOutFilePre+"_genetrees2/"+InOutFilePre+"gt2_Loci_above_25_group.txt "+ScriptFolder+" "+InOutFolder+InOutFilePre+"_genetrees2/ "+InOutFilePre+"gt2_ "+AlFileOutPost+"_al.fa fasta fasta tree "+InOutFolder+InOutFilePre+"_prunedgenetrees2 "+InOutFilePre+"pgts2_ separate 10000 10000 "+InOutFolder+InOutFilePre+"_genetrees2/"+InOutFilePre+"gt2_non_lower_outliers.txt 2 "+Mode+" "+InOutFolder+InOutFilePre+"pgt1_spptree.tre >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	Line += ScriptFolder+"tal_combiner.py "+InOutFolder+InOutFilePre+"_genetrees2/"+InOutFilePre+"gt2_Loci_above_25_group.txt "+ScriptFolder+" "+InOutFolder+InOutFilePre+"_genetrees2/ "+InOutFilePre+"gt2_ "+AlFileOutPost+"_al.fa fasta fasta "+OGName+" "+InOutFolder+InOutFilePre+"_prunedgenetrees2 "+InOutFilePre+"pgtc2_ combined 10000 10000 "+InOutFolder+InOutFilePre+"_genetrees2/"+InOutFilePre+"gt2_non_lower_outliers.txt 2 "+Mode+" >> "+InOutFolder+Date+"_"+InOutFilePre+".log\n"
	OutScript.append(Line)
	if Mode == "Parallel":
		Line = "chmod u+x "+InOutFolder+InOutFilePre+"_prunedgenetrees2/"+InOutFilePre+"pgts2_analysis_script.sh\n"
		Line += InOutFolder+InOutFilePre+"_prunedgenetrees2/"+InOutFilePre+"pgts2_analysis_script.sh\n"
		Line += "chmod u+x "+InOutFolder+InOutFilePre+"_prunedgenetrees2/"+InOutFilePre+"pgtc2_analysis_script.sh\n"
		Line += InOutFolder+InOutFilePre+"_prunedgenetrees2/"+InOutFilePre+"pgtc2_analysis_script.sh\n"
		OutScript.append(Line)
	OutFileWriting(OutFileName, OutScript)
