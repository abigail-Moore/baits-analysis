#! /usr/bin/Rscript

#modified from tlorin, StackOverflow

rm(list=ls())

#We will need ape to do phylogenetically-related calculations.
library(ape)
library(gplots)
library(RColorBrewer)

#example:
#heatmap_making.r /home/abby/transcriptomes/combined_Ln1/bestseqs_20150715/combined_alignments/RAxML_bipartitions.Ln1_comb /home/abby/transcriptomes/combined_Ln1/phot2/phot2_Loci_Obtained_Inds_reordered.txt /home/abby/transcriptomes/sandbox/test2.pdf
#heatmap_making_10groups.r TreeFileName InfoFileName OutFileName

args <- commandArgs(TRUE)
TreeFileName <- args[1]
InfoFileName <- args[2]
OutFileName <- args[3]

#TreeFileName <- '~/transcriptomes/combined2_paper1_tree/c2p1pgts2_raxmlbs/c2p1pgts2_g5_astral_rr.tre'
#InfoFileName <- '~/transcriptomes/combined2_paper1_tree/c2p1_genetrees2_text/c2p1gt2_seqs_per_ind_locus_tips_gfcombined.txt'
#OutFileName <- '~/transcriptomes/combined2_paper1_tree/c2p1_genetrees2_text/c2p1gt2_heatmap_allloci_astral_g5.pdf'


#read the tree
BackboneTree <- read.tree(TreeFileName)
is.rooted(BackboneTree)
#make it ultrametric
UMTree <- compute.brlen(BackboneTree, method="Grafen")
is.rooted(UMTree)
#make a list of the taxa to use to reorder the data matrix
taxonList <- UMTree$tip

############ get the number of taxa
#taxonList
nr1 <- length(taxonList)  ## finds the number of taxa
nr1
#########################

#turn the phylo tree to a dendrogram object
hc <- as.hclust(UMTree) #Compulsory step as as.dendrogram doesn't have a     method for phylo objects.
DGTree <- as.dendrogram(hc)
plot(DGTree, horiz=TRUE) #check dendrogram face--This does not print if I just execute the script, but does in R.

#Now reading the contig data

df <- read.table(InfoFileName, head = TRUE) #############read in your data frame

#rewrite the table in the same order (and with the same taxa) as the tree
nr <- nrow(df) ## rows in the data
nc <- ncol(df) ## columns in the data
locusNames <- colnames(df)[2:nc]#getting the names of the loci to label the final matrix
length(locusNames)

#making a new, reorganized table
df2 <- data.frame(matrix(ncol = nc, nrow = nr1)) ##make a new empty dataframe of the right size
dim(df2)
for (i in 1:nr1){    ### loop through all the taxa
taxStatus <- "not found"
for (j in 1:nr){       ## loop through the data
if (taxonList[i] == df[j,1]) ### find the row in the data equal to the ith taxon
{ 
taxStatus <- "found"
df2[i,] = df[j,]  ## put the data in he proper order, i.e., the same order as the taxon
}
}
if (taxStatus == "not found")
{
print ("ERROR!!! This taxon from tree not found in table!!!!!!")
print(taxonList[i])
}
}
dim(df2)
mat <- as.matrix(df2[,2:nc])   ### convert the data frame to a matrix
dim(mat)             #### check the dimension of the matrix
rownames(mat) <- taxonList
colnames(mat) <- locusNames


#plot the heatmaps
my_colors = c("white", "lightblue", "cornflowerblue", "darkgreen", "olivedrab", "yellow", "orange", "tomato", "darkred", "black")
#my_colors = c("white", "cornflowerblue", "darkgreen", "yellow", "orange", "darkred", "black")

pdf(OutFileName)

#plot(DGTree, horiz=TRUE) #This is just so I can plot the figure multiple times without getting erros.
par(xpd=FALSE)
heatmap.2(mat, Rowv=DGTree, Colv=NA, dendrogram='row', col = my_colors,
          sepwidth=c(0.001,0.002),sepcolor="black",colsep=1:ncol(mat)-1,rowsep=1:nrow(mat)-1,
          breaks = c(-1,0,1,2,3,4,6,8,10,15,2000),
          main = "Loci per Gene Family",
          key=FALSE,trace="none",
          cexRow=0.5,cexCol=0.4,srtCol=45,
          margins=c(5,10))
par(lend=2)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("none", "1", "2", "3", "4", "5-6", "7-8", "9-10", "11-15", ">15"),
       fill = my_colors,  # color key
       border = "black",
       ncol = 4
)

dev.off()
