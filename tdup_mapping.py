#! /usr/bin/env python

#tdup_mapping.py v. 1.0 19 April 2016 Abby Moore
#This script takes a species tree and the list of positions of the duplications, produced by
#tnotung_homolog_parsing.py, and maps those on the tree, producing a tree with the
#number of duplications per node.
#If the position of a duplication is uncertain, then it is divided equally between
#the possible nodes (so if a duplication could be at 2 nodes, each node gets half
#a duplication

#dup_pos_dict.txt (produced by tnotung_homolog_parsing.py)
'''
First line is header, tab-delimitted
Locus [0]: mdh
Duplication_Name [1]: Y28
Duplicated_Individuals [2] (comma separated): Calandrinia_crispisepala_21
Sister_Group [3] (comma separated): Calandrinia_translucens_99,Calandrinia_translucens_4,Calandrinia_lehmannii_35,Calandrinia_kalenninsis_8
'''

import sys
import re
from collections import defaultdict, OrderedDict
import dendropy

'''
tdup_mapping.py ~/transcriptomes/combined_Ln12_tr/Ln12s_prunedgenetrees/RAxML_bipartitions.Ln12spgt_combined_64inds_185seqs_rooted ~/transcriptomes/combined_Ln12_tr/tomovespades/Ln12sgt_dup_pos_dict_head.txt
tdup_mapping.py ~/transcriptomes/general/spp_tree_Ln12_ts_comb_nn ~/transcriptomes/combined_Ln12_tr/tomovespades/Ln12sgt_dup_pos_dict.txt ~/transcriptomes/combined_Ln12_tr/tomovespades/Ln12sgt_tree_dups_mapped
tdup_mapping.py SppTreeFN DPDFN OutFileName DecPlaces[none, two]
'''

Usage = '''
tdup_mapping.py version 1.0
This is a script to map the locations of duplications on a tree.  Or that is 
what it will theoretically be, if it works....
tdup_mapping.py
[species tree file name]
[dup_pos_dict.txt file produced by tnotung_homolog_parsing.py]
[name for the output tree file with number of duplications mapped as node 
labels]
[number of decimal places to be shown for duplication node labels, "none" for 
none (what FigTree needs) or "two" for 2]
'''

DecPlacesList = ["none", "two"]

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 5:
	sys.exit("ERROR!  This script requires 4 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
SppTreeFN = sys.argv[1]
DPDFN = sys.argv[2]
OutFileName = sys.argv[3]
DecPlaces = sys.argv[4]
if DecPlaces not in DecPlacesList:
	sys.exit ("ERROR!  The number of decimal places can only be %s, and you wrote %s.\n%s" % (", ".join(DecPlacesList), DecPlaces, Usage))

Verbose = False
ShowTrees = False

Verbose = True
ShowTrees = True

#################################################################################

#NodeDistFinder looks for the nodes that are between two nodes when the second is directly
#ancestral to the first.
#original
def NodeDistFinder(TempTree, Node1, Node2):
	if Node1 == Node2:
		ListTemp = [ ]
		NN = 0
	else:
		ListTemp = [ Node1 ]
		NN = 1
		for ANode in Node1.ancestor_iter():
			if ANode == Node2:
				break
			else:
				ListTemp.append(ANode)
				NN += 1
	return (ListTemp, NN)
	#These are the various (ListIntNodes, NumIntNodes).

#IndListfromTree makes a list of the sequences present in a tree
#This is no longer being used, but I will leave it for now, since this is the only script it is in.
#from tcombparse_to_trees.py
def IndListfromTree(FileName):
	InFile = open(FileName, 'rU')
	for Line in InFile:
		TempTree = dendropy.Tree.get_from_string(Line.strip('\n').strip('\r'), schema = 'newick', preserve_underscores=True)
	InFile.close()
	ListTemp = [Node.taxon.label for Node in TempTree.leaf_nodes()]
	return ListTemp
	#This is IndList.

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
	return

#################################################################################

#read the species tree
SppTree = dendropy.Tree.get(path=SppTreeFN, schema='newick', preserve_underscores=True)
#if ShowTrees == True: print(SppTree.as_ascii_plot(show_internal_node_labels=True))
#read the species list from the tree
IndList = IndListfromTree(SppTreeFN)

#read the dup_pos_dict
DPDO = defaultdict(dict)
DPDFile = open(DPDFN, 'rU')
for Line in DPDFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	if Line[0] != "Locus":
		DupName = Line[0]+"_"+Line[1]
		InGroupList = Line[2].split(',')
		SisterGroupList = Line[3].split(',')
		DPDO[DupName]['InGroup'] = InGroupList
		DPDO[DupName]['SisterGroup'] = SisterGroupList
DPDFile.close()

#making a revised DPD with only the individuals found in the tree
DPD = defaultdict(dict)
for DupName in DPDO:
	InGroupList = list(set(DPDO[DupName]['InGroup']).intersection(IndList))
	InGroupOut = list(set(DPDO[DupName]['InGroup']).difference(IndList))
	if len(InGroupOut) != 0:
		print('The following ingroup sequences were not in the species tree: %s.\n' % (', '.join(InGroupOut)))
	SisterGroupList = list(set(DPDO[DupName]['SisterGroup']).intersection(IndList))
	if (len(InGroupList) != 0) and (len(SisterGroupList) != 0):
		DPD[DupName]['InGroup'] = InGroupList
		DPD[DupName]['SisterGroup'] = SisterGroupList

#for each duplication:
for DupName in DPD:
	#find the node subtending all duplicated individuals
	DupNode = SppTree.mrca(taxon_labels = DPD[DupName]['InGroup'])
	#find the node subtending all sister individuals
	SisNode = SppTree.mrca(taxon_labels = DPD[DupName]['SisterGroup'])
	PostDupNode = SppTree.mrca(taxon_labels = list(set(DPD[DupName]['InGroup']+DPD[DupName]['SisterGroup'])))
	(ListIntNodes, NumIntNodes) = NodeDistFinder(SppTree, DupNode, PostDupNode)
	if Verbose == True: print("There are %d nodes between the duplication and its mrca with its sister group for duplication %s.\n" % (NumIntNodes, DupName))
	#There are various possibilities for how the placement of the duplication:
	#if the DupNode is separated by one from the PostDupNode, the duplication is mapped to the DupNode.
	if NumIntNodes == 1:
		if Verbose == True: print("Duplication %s is resolved to one node, as the sister group is at the sister node to the duplication!\n" % (DupName))
		try:
			DupNode.label += 1.0
		except TypeError:
			DupNode.label = 1.0
	#if the DupNode and the PostDupNode are the same, then there were two sequential duplications at the same node, and the duplication is also mapped to the
	#DupNode.
	elif NumIntNodes == 0:
		if Verbose == True: print("Duplication %s is resolved to one node, as the sister node is the same as the duplicated node.\n" % (DupName))
		try:
			DupNode.label += 1.0
		except TypeError:
			DupNode.label = 1.0
	#if the DupNode is separated by multiple nodes from the PostDupNode (because the intermediate individuals were not present in that gene tree),
	#the duplication is spread evenly over the intervening nodes.
	elif NumIntNodes > 1:
		if Verbose == True:
			print("There are %d intermediate nodes between duplication %s and its sister group.\n" % (NumIntNodes, DupName))
			print("Each of these nodes will be given %f of a duplication.\n" % (1.0/NumIntNodes))
		for IntNode in ListIntNodes:
			try:
				IntNode.label += 1.0/NumIntNodes
			except TypeError:
				IntNode.label = 1.0/NumIntNodes

for Node in SppTree.preorder_node_iter():
	if Node.label == None:
		Node.label = "0"
	else:
		#This is where you decide if you want whole numbers for the node labels and, if not, how many decimal places you want the labels to have.
		#However, figtree has a problem reading node labels with decimal places.
		#node labels without decimal places
		if DecPlaces == "none":
			Node.label = str(int(round(float(Node.label))))
		#node labels with decimal places:
		elif DecPlaces == "two":
			if ((Node.label % 1) == 0):
				Node.label = str(int(Node.label))
			else:
				Node.label = ("%.2f" % (Node.label))
		#print Node.label
		
if ShowTrees == True: print(SppTree.as_ascii_plot(show_internal_node_labels=True))

SppTree.write(path=OutFileName, schema="newick", suppress_leaf_node_labels = False, suppress_internal_node_labels = False, unquoted_underscores=True, node_label_element_separator = ",")

