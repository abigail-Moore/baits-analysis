#! /usr/bin/Rscript

#tfigmaking.r version 1.0 16 Feb. 2017 Abby Moore
#This script makes various figures for bait coverage from the python output.

rm(list=ls())

args <- commandArgs(TRUE)
InFolder <- args[1]
InFilePre <- args[2]
OutFolder <- args[3]

#InFolder = "/home/abby/transcriptomes/combined2_paper1_tree/c2p1_genetrees2_text/"
#InFilePre = "c2p1gt2_"
#OutFolder = "/home/abby/transcriptomes/combined2_paper1_tree/c2p1_genetrees2_text/"

#reading the data
InFileName <- paste(InFolder, InFilePre, "length_vs_numinds_all.csv", sep = "")
InfoTable <- read.csv(InFileName, header = TRUE, row.names = 1)
str(InfoTable)
#row.names(InfoTable)

#sequence length versus number of individuals
OutFileName <- paste(OutFolder, InFilePre,'NumInds_vs_MeanSeqLen.pdf', sep = "")
pdf(OutFileName)
plot (x=InfoTable$num_inds, y = InfoTable$mean_length, xlab = "Number of Individuals", ylab = "Mean Sequence Length")
abline(lm(InfoTable$mean_length~InfoTable$num_inds))
summary(lm(InfoTable$mean_length~InfoTable$num_inds))
dev.off()

###number of individuals alone
#sorting
NumIndsSorted <- sort(InfoTable$num_inds, decreasing = TRUE)
#making a vector that is the same length as the number of loci, so we can plot them
NumLoci <- length(NumIndsSorted)
LocusNumList <- c(1:NumLoci)
OutFileName <- paste(OutFolder, InFilePre,'NumInds_per_Locus.pdf', sep = "")
pdf(OutFileName)
plot(x = LocusNumList, y = NumIndsSorted, type = "l", xlab = "Locus Rank (by individual)", ylab = "Number of Individuals", main = "Individuals per Locus")
dev.off()

##locus length alone
MeanLengthSorted <- sort(InfoTable$mean_length, decreasing = TRUE)
OutFileName <- paste(OutFolder, InFilePre,'MeanLocusLength.pdf', sep = "")
pdf(OutFileName)
plot(x = LocusNumList, y = MeanLengthSorted, type = "l", xlab = "Locus Rank (by length)", ylab = "Mean Length", main = "Mean Locus Length")
dev.off()

###attempting to count the number of duplications per gene family
#InFileName2 <- paste(InFolder, InFilePre, 'dup_pos_dict.txt', sep = "")
#InfoTable2 <- read.table(InFileName2, header = TRUE, sep = "\t")
#str(InfoTable2)

#counting the number of duplications per gene family
#DupsPerGF <- tapply(InfoTable2$Locus,InfoTable2$Locus,length)
#write.table(DupsPerGF, file = "~/transcriptomes/combined2_paper1_tree/c2p1_genetrees2_txt/c2p1gt2_dups.txt", sep = "\t")

#sorting
#sortedDpGF <- sort(DupsPerGF, decreasing = TRUE)
#adding one, because the number of loci is the number of duplications plus one.
#LociPerGF <- sortedDpGF + 1
#str(sortedDpGF)
#str(LociPerGF)
#NumGFs <- length(DupsPerGF)
#GFNumList <- c(1:NumGFs)
#OutFileName <- paste(OutFolder, InFilePre,'Loci_per_GF.pdf', sep = "")
#pdf(OutFileName)
#plot(x = GFNumList, y = LociPerGF, type = "l", xlab = "Gene Family Rank (by number of loci)", ylab = "Number of Loci", main = "Number of Duplications per Gene Family")
#dev.off()

#reading in the new table
InFileName3 <- paste(InFolder, InFilePre, 'dups.csv', sep = "")
InfoTable3 <- read.csv (InFileName3, header = TRUE)
str(InfoTable3)
#sorting
sortedLocipGF <- sort(InfoTable3$Num_Copies)
NumGFs <- length(sortedLocipGF)
GFNumList <- c(1:NumGFs)
OutFileName <- paste(OutFolder, InFilePre,'Loci_per_GF.pdf', sep = "")
pdf(OutFileName)
plot(x = GFNumList, y = sortedLocipGF, type = "l", xlab = "Gene Families", ylab = "Number of Duplications", main = "Number of Duplications per Gene Family")
dev.off()


#number of individuals versus number of groups per locus

InFileName4 <- paste(InFolder, InFilePre, 'Inds_per_Paralog.csv', sep = "")
InfoTable4 <- read.csv (InFileName4, header = TRUE)
str(InfoTable4)
#just families per locus for now
sortedNumGroups <- sort(InfoTable4$Number_of_Groups_with_Sequences)
OutFileName <- paste(OutFolder, InFilePre, 'Groups_per_Locus.pdf', sep = "")
pdf(OutFileName)
plot(x = LocusNumList, y = sortedNumGroups, type = "l", xlab = "Locus", ylab = "Number of Groups", main = "Number of Groups per Locus")
dev.off()

#then numinds versus numfams
OutFileName <- paste(OutFolder, InFilePre, 'Inds_versus_Groups_per_Locus.pdf', sep = "")
pdf(OutFileName)
plot(x = InfoTable4$Total_Sequences_in_All_Groups, y = InfoTable4$Number_of_Groups_with_Sequences, xlab = "Number of Individuals with that Locus", ylab = "Number of Groups with that Locus")
abline(lm(InfoTable4$Number_of_Groups_with_Sequences~InfoTable4$Total_Sequences_in_All_Groups))
summary(lm(InfoTable4$Number_of_Groups_with_Sequences~InfoTable4$Total_Sequences_in_All_Groups))
dev.off()

#We concluded that we did not need a graph of the number of individuals versus the length
#of the alignment, so I tfig_making_prep.py does not make inds_seqlen.txt, so this does not work.
#InFileName5 <- paste(InFolder, InFilePre, 'inds_seqlen.txt', sep = "")
#InfoTable5 <- read.table(InFileName5, header = TRUE, sep = "\t")
#str(InfoTable5)
#NumLociUsed <- length(InfoTable5$Locus)
#NumLociUsed
#NumLociUsedList <- c(1:NumLociUsed)

#OutFileName <- paste(OutFolder, InFilePre,'NumInds_vs_AlLen.pdf', sep = "")
#pdf(OutFileName)
#plot (x=InfoTable5$NumInds, y = InfoTable5$SeqLen, xlab = "Number of Individuals", ylab = "Alignment Length")
#abline(lm(InfoTable5$SeqLen~InfoTable5$NumInds))
#summary(lm(InfoTable5$SeqLen~InfoTable5$NumInds))
#dev.off()

#sortedAlLength <- sort(InfoTable5$SeqLen)
#OutFileName <- paste(OutFolder, InFilePre, 'Alignment_Length.pdf', sep = "")
#pdf(OutFileName)
#plot(x = NumLociUsedList, y = sortedAlLength, type = "l", xlab = "Locus", ylab = "Alignment Length", main = "Length of Sequence Alignment per Locus")
#dev.off()
